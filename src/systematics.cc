#include <systematics.h>
#include <types.h>
#include <utilities.h>
#include <commandlinecfg.h>

#include "yaml-cpp/yaml.h"

#include <Math/PdfFunc.h>
#include <TH1.h>

#include <algorithm>
#include <iostream>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace plotIt {
    SystematicSet::SystematicSet(Systematic& parent):
        parent(&parent) {

    }

    void SystematicSet::update() {
        parent->apply(*this);
    }

    void SystematicSet::scale(float factor) {
        if (nominal_shape) {
            static_cast<TH1*>(nominal_shape.get())->Scale(factor);
        }

        if (up_shape) {
            static_cast<TH1*>(up_shape.get())->Scale(factor);
        }

        if (down_shape) {
            static_cast<TH1*>(down_shape.get())->Scale(factor);
        }
    }

    void SystematicSet::scale(float factor, std::string scale_option_syst) {
        if (nominal_shape) {
            static_cast<TH1*>(nominal_shape.get())->Scale(factor, scale_option_syst.c_str());
        }

        if (up_shape) {
            static_cast<TH1*>(up_shape.get())->Scale(factor, scale_option_syst.c_str());
        }

        if (down_shape) {
            static_cast<TH1*>(down_shape.get())->Scale(factor, scale_option_syst.c_str());
        }
    }

    void SystematicSet::rebin(size_t factor) {
        if (nominal_shape) {
            static_cast<TH1*>(nominal_shape.get())->Rebin(factor);
        }

        if (up_shape) {
            static_cast<TH1*>(up_shape.get())->Rebin(factor);
        }

        if (down_shape) {
            static_cast<TH1*>(down_shape.get())->Rebin(factor);
        }
    }

    std::string SystematicSet::name() const {
        return parent->name;
    }

    std::string SystematicSet::prettyName() const {
        return parent->pretty_name;
    }

    SystematicSet Systematic::newSet(TObject* nominal, File& file, const Plot& plot) {
        SystematicSet s = SystematicSet(*this);
        s.true_nominal_shape.reset(nominal->Clone());
        s.true_up_shape.reset(nominal->Clone());
        s.true_down_shape.reset(nominal->Clone());

        return s;
    }

    void Systematic::apply(SystematicSet& systs) {
        systs.nominal_shape.reset(systs.true_nominal_shape->Clone());
        systs.up_shape.reset(systs.true_up_shape->Clone());
        systs.down_shape.reset(systs.true_down_shape->Clone());
    }

    ConstantSystematic::ConstantSystematic(const YAML::Node& node) {
        if (node.IsScalar()) {
            value = node.as<float>();
            return;
        }

        value = node["value"].as<float>();
    }

    void ConstantSystematic::apply(SystematicSet& systs) {
        Systematic::apply(systs);

        TH1* up = static_cast<TH1*>(systs.up_shape.get());
        TH1* down = static_cast<TH1*>(systs.down_shape.get());

        up->Scale(value);
        down->Scale(2 - value);
    }

    LogNormalSystematic::LogNormalSystematic(const YAML::Node& node) {
        if (node.IsScalar()) {
            prior = node.as<float>();
            eval();

            return;
        }

        prior = node["prior"].as<float>();
        if (node["post-fit"])
            postfit = node["post-fit"].as<float>();

        if (node["post-fit-error"]) {
            postfit_error_up = postfit_error_down = node["post-fit-error"].as<float>();
        }

        if (node["post-fit-error-up"]) {
            postfit_error_up = node["post-fit-error-up"].as<float>();
        }

        if (node["post-fit-error-down"]) {
            postfit_error_down = node["post-fit-error-down"].as<float>();
        }

        eval();
    }

    void LogNormalSystematic::apply(SystematicSet& systs) {
        Systematic::apply(systs);

        TH1* up = static_cast<TH1*>(systs.up_shape.get());
        TH1* down = static_cast<TH1*>(systs.down_shape.get());

        up->Scale(value_up);
        down->Scale(value_down);
    }

    void LogNormalSystematic::eval() {
        value = exp(postfit * log(prior));
        value_up = exp((postfit + postfit_error_up) * log(prior));
        value_down = exp((postfit - postfit_error_down) * log(prior));
    }

    ShapeSystematic::ShapeSystematic(const YAML::Node& node) {

        if(node["ext-sum-weight-up"])
            ext_sum_weight_up = node["ext-sum-weight-up"].as<float>();
        if(node["ext-sum-weight-down"])
            ext_sum_weight_down = node["ext-sum-weight-down"].as<float>();

    }

    SystematicSet ShapeSystematic::newSet(TObject* nominal, File& file, const Plot& plot) {

        auto result = Systematic::newSet(nominal, file, plot);

        // We need to find the up and down shape
        // Two possibilities:
        //   - we look for two objects named <nominal>__<systematic>[up|down] in the same file
        //   - we look for two objects named <nominal> in the file <nominal>__<systematic>[up|down].root

        std::array<Variation, 2> variations = {UP, DOWN};
        std::map<Variation, std::shared_ptr<TObject>*> links = {{UP, &result.true_up_shape}, {DOWN, &result.true_down_shape}};

        auto formatSystematicsName = [this](Variation variation) {
            static std::map<Variation, std::string> names = {{UP, "up"}, {DOWN, "down"}};
            return "__" + this->name + names[variation];
        };

        auto formatSystematicsName2 = [this](Variation variation) {
            static std::map<Variation, std::string> names = {{UP, "Up"}, {DOWN, "Down"}};
            return "__" + this->name + names[variation];
        };

        for (const auto& variation: variations) {
            std::string object_postfix = formatSystematicsName(variation);

            TObject* object = NULL;

            if (!CommandLineCfg::get().desytop) {

                std::string object_name = applyRenaming(file.renaming_ops, plot.name) + object_postfix;
                object = file.handle->Get(object_name.c_str());

                if (!object) {
                    std::string object_postfix2 = formatSystematicsName2(variation);
                    std::string object_name2 = applyRenaming(file.renaming_ops, plot.name) + object_postfix2;
                    object = file.handle->Get(object_name2.c_str());
                }

                if (object) {
                    links[variation]->reset(object->Clone());
                    continue;
                }
            }

            auto nominal_path = fs::path(file.path);
            auto syst_path = nominal_path.parent_path();
            syst_path /= nominal_path.stem();
            syst_path += object_postfix;
            syst_path += ".root";

            if (fs::exists(syst_path)) {
                std::shared_ptr<TFile>& f = file.friend_handles[syst_path.native()];
                if (! f)
                    f.reset(TFile::Open(syst_path.native().c_str()));

                if (ext_sum_weight_up > 1.1 and ext_sum_weight_down > 1.1){
                    if(variation == UP) {
                        TH1F* tmp = (TH1F*) f->Get(plot.name.c_str());
                        tmp->Scale(file.generated_events / ext_sum_weight_up);
                        object = tmp;
                    }
                    else if (variation == DOWN) {
                        TH1F* tmp = (TH1F*) f->Get(plot.name.c_str());
                        tmp->Scale(file.generated_events / ext_sum_weight_down);
                        object = tmp;
                    }
                }

                else object = f->Get(plot.name.c_str());

                if (object) {

                    // Rebin by array
                    std::shared_ptr<TObject> cloned_obj(object->Clone());
                    if (std::size(plot.rebin_arr) > 0) {
                        auto h = std::static_pointer_cast<TH1>(cloned_obj).get();
                        auto size_arr = std::size(plot.rebin_arr);
                        double xbins[size_arr];
                        for (size_t i=0; i<size_arr; i++) {
                            xbins[i] = plot.rebin_arr.at(i);
                        }
                        auto hnew = h->Rebin(int(size_arr)-1, plot.name.c_str(), xbins);
                        cloned_obj.reset(hnew);
                        object = cloned_obj.get();
                    }

                    links[variation]->reset(object->Clone());
                }
            }
        }

        return result;
    }


    std::shared_ptr<Systematic> SystematicFactory::create(const std::string& name, const std::string& type, const YAML::Node& node) {

        std::string lower_type = type;
        std::transform(lower_type.begin(), lower_type.end(), lower_type.begin(), ::tolower);

        std::shared_ptr<Systematic> result;
        if (type == "const") {
            result = std::make_shared<ConstantSystematic>(node);
        } else if (type == "lognormal" || type == "ln") {
            result = std::make_shared<LogNormalSystematic>(node);
        } else if (type == "shape") {
            result = std::make_shared<ShapeSystematic>(node);
        } else if (type == "siglike") {
            result = std::make_shared<ShapeSystematic>(node);
        }

        if (result) {
            result->name = name;
            result->pretty_name = name;

            std::string on = ".*";
            if (node.IsMap()) {
                if (node["pretty-name"])
                    result->pretty_name = node["pretty-name"].as<std::string>();

                if (node["on"])
                    on = node["on"].as<std::string>();
            }

            result->on = std::regex(on);

            return result;
        }

        throw std::invalid_argument("Unknown systematic type: " + type);
    }
};
