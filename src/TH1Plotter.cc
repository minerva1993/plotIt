#include <TH1Plotter.h>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <TLine.h>
#include <TObject.h>
#include <TPave.h>
#include <TVirtualFitter.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TLegend.h>
#include <TROOT.h>

#include <commandlinecfg.h>
#include <pool.h>
#include <utilities.h>

namespace plotIt {

    /*!
     * Compute the ratio between two histograms, taking into account asymetric error bars
     */
    std::shared_ptr<TGraphAsymmErrors> getRatio(TH1* a, TH1* b) {
        std::shared_ptr<TGraphAsymmErrors> g(new TGraphAsymmErrors(a));

        size_t npoint = 0;
        for (size_t i = 1; i <= (size_t) a->GetNbinsX(); i++) {
            float b1 = a->GetBinContent(i);
            float b2 = b->GetBinContent(i);

            if ((b1 == 0) || (b2 == 0))
                continue;

            float ratio = b1 / b2;

            float b1sq = b1 * b1;
            float b2sq = b2 * b2;

            float e1sq_up = a->GetBinErrorUp(i) * a->GetBinErrorUp(i);
            float e2sq_up = b->GetBinErrorUp(i) * b->GetBinErrorUp(i);

            float e1sq_low = a->GetBinErrorLow(i) * a->GetBinErrorLow(i);
            float e2sq_low = b->GetBinErrorLow(i) * b->GetBinErrorLow(i);

            float error_up = sqrt((e1sq_up * b2sq + e2sq_up * b1sq) / (b2sq * b2sq));
            float error_low = sqrt((e1sq_low * b2sq + e2sq_low * b1sq) / (b2sq * b2sq));

            //Set the point center and its errors
            g->SetPoint(npoint, a->GetBinCenter(i), ratio);
            //g->SetPointError(npoint, a->GetBinCenter(i) - a->GetBinLowEdge(i),
                    //a->GetBinLowEdge(i) - a->GetBinCenter(i) + a->GetBinWidth(i),
                    //error_low, error_up);
            g->SetPointError(npoint, 0, 0, error_low, error_up);
            npoint++;
        }

        g->Set(npoint);

        return g;
    }


    /*Not to propagate mc stat unc to data, use this function*/
    std::shared_ptr<TGraphAsymmErrors> getRatio2(TH1* a, TH1* b) {
        std::shared_ptr<TGraphAsymmErrors> g(new TGraphAsymmErrors(a));

        size_t npoint = 0;
        for (size_t i = 1; i <= (size_t) a->GetNbinsX(); i++) {
            float b1 = a->GetBinContent(i);
            float b2 = b->GetBinContent(i);

            if ((b1 == 0) || (b2 == 0))
                continue;

            float ratio = b1 / b2;

            float error_up = std::abs(a->GetBinErrorUp(i) / b2);
            float error_low = std::abs(a->GetBinErrorLow(i) / b2);

            //Set the point center and its errors
            g->SetPoint(npoint, a->GetBinCenter(i), ratio);
            g->SetPointError(npoint, 0, 0, error_low, error_up);
            npoint++;
        }

        g->Set(npoint);

        return g;
    }

  bool TH1Plotter::supports(TObject& object) {
    return object.InheritsFrom("TH1");
  }

  TH1Plotter::Stacks TH1Plotter::buildStacks(bool sortByYields) {
      std::set<int64_t> indices;

      for (auto& file: m_plotIt.getFiles()) {
          if (file.type == MC)
              indices.emplace(file.stack_index);
      }

      Stacks stacks;
      for (auto index: indices) {
          auto stack = buildStack(index, sortByYields);
          if (stack.stack)
              stacks.push_back(std::make_pair(index, stack));
      }

      return stacks;
  }

  TH1Plotter::Stack TH1Plotter::buildStack(int64_t index, bool sortByYields) {

      std::shared_ptr<THStack> stack;
      std::shared_ptr<TH1> histo_merged;

      std::string stack_name = "mc_stack_" + std::to_string(index);

      // First pass. Merged all member of a group into a single
      // histogram.
      // Key is group name, value is group histogram
      std::vector<std::pair<std::string, std::shared_ptr<TH1>>> group_histograms;
      for ( auto& file: m_plotIt.getFiles([this,index] ( const File& f ) {
            return ( f.type == MC ) && ( f.stack_index == index )
                && ( ! f.legend_group.empty() ) && ( dynamic_cast<TH1*>(f.object)->GetEntries() != 0 );
            } ) ) {
          TH1* nominal = dynamic_cast<TH1*>(file.object);

          // Remove small mc bins...always very small number is filled.
          if (CommandLineCfg::get().desytop) {
              for (int i=1; i <= nominal->GetNbinsX(); i++) {
                  float tmpbin = nominal->GetBinContent(i);
                  //in code, 0.0001 is filled, then samples are merged.
                  if (tmpbin < 0.005) {
                      nominal->SetBinContent(i, 0.);
                      nominal->SetBinError(i, 0.);
                  }
              }
          }

          auto it = std::find_if(group_histograms.begin(), group_histograms.end(), [&file](const std::pair<std::string, std::shared_ptr<TH1>>& item) {
                  return item.first == file.legend_group;
                  });

          if (it == group_histograms.end()) {
              std::string name = "group_histo_" + file.legend_group + "_" + stack_name;
              std::shared_ptr<TH1> h(dynamic_cast<TH1*>(nominal->Clone(name.c_str())));
              h->SetDirectory(nullptr);
              group_histograms.push_back(std::make_pair(file.legend_group, h));
          } else {
              it->second->Add(nominal);
          }
      }

      std::vector<std::tuple<TH1*, std::string>> histograms_in_stack;

      for ( auto& file: m_plotIt.getFiles([this,index] ( const File& f ) {
            return ( f.type == MC ) && ( f.stack_index == index )
                && ( ! ( ( dynamic_cast<TH1*>(f.object)->GetEntries() == 0 ) && (f.legend_group.empty()) ) );
            } ) ) {

          TH1* nominal = dynamic_cast<TH1*>(file.object);
          if (!stack) {
              stack = std::make_shared<THStack>(stack_name.c_str(), stack_name.c_str());
              TemporaryPool::get().add(stack);
          }

          // Try to find if this file is a member of a group
          // If it is, then use the merged histogram built above, but only for
          // the first file.
          auto it = std::find_if(group_histograms.begin(), group_histograms.end(), [&file](const std::pair<std::string, std::shared_ptr<TH1>>& item) {
                  return item.first == file.legend_group;
                  });

          if (!file.legend_group.empty() && it == group_histograms.end()) {
              // The group histogram has already been added to the stack, so
              // skip to the next one
              continue;
          } else if (it != group_histograms.end()) {
              auto n = it->second;
              TemporaryPool::get().add(n);

              nominal = n.get();

              // Since we will add this group histogram to the stack
              // remove it from the pool to avoid double addition
              group_histograms.erase(it);
          }

          histograms_in_stack.push_back(std::make_tuple(nominal, m_plotIt.getPlotStyle(file)->drawing_options));
      }

      // Sort histograms by yields
      if (sortByYields) {
          std::sort(histograms_in_stack.begin(), histograms_in_stack.end(),
                  [](const std::tuple<TH1*, std::string>& a, const std::tuple<TH1*, std::string>& b) {
                    return std::get<0>(a)->Integral() < std::get<0>(b)->Integral();
                  });
      }

      for (const auto& t: histograms_in_stack) {
          TH1* nominal = std::get<0>(t);
          stack->Add(nominal, std::get<1>(t).c_str());

          if (histo_merged) {
              histo_merged->Add(nominal);
          } else {
              std::string name = "mc_stat_only_" + stack_name;
              histo_merged.reset( dynamic_cast<TH1*>(nominal->Clone(name.c_str())) );
              histo_merged->SetDirectory(nullptr);
          }
      }

      // Ensure there's MC events
      if (stack && !stack->GetHists()) {
          histo_merged.reset();
          stack.reset();
      }

      Stack s {stack, histo_merged};

      return s;
  }

  void TH1Plotter::computeSystematics(int64_t index, Stack& stack, Summary& summary) {

      // Key is systematics name, value is the combined systematics value for each bin
      std::map<std::string, std::vector<float>> combined_systematics_map;
      std::map<std::string, std::vector<float>> combined_systematics_map_up;
      std::map<std::string, std::vector<float>> combined_systematics_map_dn;

      for ( auto& file: m_plotIt.getFiles([this,index] ( const File& f ) {
            return ( f.type != DATA ) && ( ! f.systematics->empty() )
                && ( ( f.type != MC ) || ( f.stack_index == index ) ) ;
            } ) ) {

          for (auto& syst: *file.systematics) {

              std::vector<float>* combined_systematics;
              std::vector<float>* combined_systematics_up;
              std::vector<float>* combined_systematics_dn;
              auto map_it = combined_systematics_map.find(syst.name());
              auto map_it_up = combined_systematics_map_up.find(syst.name());
              auto map_it_dn = combined_systematics_map_dn.find(syst.name());

              if (map_it == combined_systematics_map.end()) {
                  combined_systematics = &combined_systematics_map[syst.name()];
                  combined_systematics->resize(stack.syst_only->GetNbinsX(), 0.);
              } else {
                  combined_systematics = &map_it->second;
              }
              if (map_it_up == combined_systematics_map_up.end()) {
                  combined_systematics_up = &combined_systematics_map_up[syst.name()];
                  combined_systematics_up->resize(stack.syst_only->GetNbinsX(), 0.);
              } else {
                  combined_systematics_up = &map_it_up->second;
              }
              if (map_it_dn == combined_systematics_map_dn.end()) {
                  combined_systematics_dn = &combined_systematics_map_dn[syst.name()];
                  combined_systematics_dn->resize(stack.syst_only->GetNbinsX(), 0.);
              } else {
                  combined_systematics_dn = &map_it_dn->second;
              }

              TH1* nominal_shape = static_cast<TH1*>(syst.nominal_shape.get());
              TH1* up_shape = static_cast<TH1*>(syst.up_shape.get());
              TH1* down_shape = static_cast<TH1*>(syst.down_shape.get());

              if (! nominal_shape || ! up_shape || ! down_shape)
                  continue;


              float total_syst_error = 0;
              float total_syst_error_up = 0;
              float total_syst_error_dn = 0;

              // First, calculate up/down yield variation
              // here, up is defined by the integram of (var-nom) > 0
              // If one-sided, uncs are square-summed
              float temp_total_syst_error_up = up_shape->Integral() - nominal_shape->Integral();
              float temp_total_syst_error_down = down_shape->Integral() - nominal_shape->Integral();
              if (temp_total_syst_error_up * temp_total_syst_error_down <= 0) {
                  total_syst_error_up = temp_total_syst_error_up;
                  total_syst_error_dn = temp_total_syst_error_down;
              } else {
                  total_syst_error_up = temp_total_syst_error_down;
                  total_syst_error_dn = temp_total_syst_error_up;
              }

              // Systematics in each bin are fully correlated, as they come either from
              // a global variation, or for a shape variation. The total systematics error
              // is simply for sum of all errors in each bins
              //
              // However, we consider that different systematics in the same bin are totaly
              // uncorrelated. The total systematics errors is then the quadratic sum.
              // Thus, first calculate linear sum of variation per systematic source,
              // then square sum by looping over combined_systematics_map
              //
              // In addition, we take only larger variation in the case of one-sided unc.
              for (uint32_t i = 1; i <= (uint32_t) stack.syst_only->GetNbinsX(); i++) {

                  float syst_error_up = 0.;
                  float syst_error_dn = 0.;
                  float temp_syst_error_up = up_shape->GetBinContent(i) - nominal_shape->GetBinContent(i);
                  float temp_syst_error_dn = down_shape->GetBinContent(i) - nominal_shape->GetBinContent(i);

                  // Remove very small variations by hand - not applied for now
                  //if (std::abs(temp_syst_error_up) < 1e-07) temp_syst_error_up = 1e-07;
                  //if (std::abs(temp_syst_error_dn) < 1e-07) temp_syst_error_dn = 1e-07;

                  // Normal case, double-sided
                  if (temp_syst_error_up * temp_syst_error_dn <= 0) {
                      if (temp_syst_error_up > 0 and temp_syst_error_dn <= 0) {
                          syst_error_up = temp_syst_error_up;
                          syst_error_dn = temp_syst_error_dn;
                      } else {
                          syst_error_up = temp_syst_error_dn;
                          syst_error_dn = temp_syst_error_up;
                      }
                  // One-sided.
                  } else {
                      if (temp_syst_error_up > 0) {
                          syst_error_up = std::max(temp_syst_error_up, temp_syst_error_dn);
                          syst_error_dn = 0.0;
                      } else {
                          syst_error_up = 0.0;
                          syst_error_dn = std::max(temp_syst_error_up, temp_syst_error_dn);
                      }
                  }
                  float syst_error = std::max(syst_error_up, syst_error_dn);

                  total_syst_error += syst_error;
                  //std::cout << "per bin " << syst.name() << " " << syst_error << " " << syst_error_up << " " << syst_error_dn << std::endl;

                  // Only propagate uncertainties for MC, not signal
                  if (file.type == MC) {
                      (*combined_systematics)[i - 1] += syst_error;
                      (*combined_systematics_up)[i - 1] += syst_error_up;
                      (*combined_systematics_dn)[i - 1] += syst_error_dn;
                  }
              }

              SummaryItem summary_item;
              summary_item.process_id = file.id;
              summary_item.name = syst.prettyName();
              summary_item.events_uncertainty = total_syst_error;
              summary_item.events_uncertainty_up = total_syst_error_up;
              summary_item.events_uncertainty_down = total_syst_error_dn;

              summary.addSystematics(file.type, file.id, summary_item);
          }
      }

      // Combine all systematics in one
      // Consider that all the systematics are not correlated
      for (auto& combined_systematics: combined_systematics_map) {
          for (size_t i = 1; i <= (size_t) stack.syst_only->GetNbinsX(); i++) {
              float total_error = stack.syst_only->GetBinError(i);
              stack.syst_only->SetBinError(i, std::sqrt(total_error * total_error + combined_systematics.second[i - 1] * combined_systematics.second[i - 1]));
          }
      }

      for (auto& combined_systematics_up: combined_systematics_map_up) {
          for (size_t i = 0; i < (size_t) stack.syst_only->GetNbinsX(); i++) {
              float total_error_up = stack.syst_only_asym->GetErrorYhigh(i);
              stack.syst_only_asym->SetPointEYhigh(i, std::sqrt(total_error_up * total_error_up + combined_systematics_up.second[i] * combined_systematics_up.second[i]));
          }
      }

      for (auto& combined_systematics_dn: combined_systematics_map_dn) {
          for (size_t i = 0; i < (size_t) stack.syst_only->GetNbinsX(); i++) {
              float total_error_down = stack.syst_only_asym->GetErrorYlow(i);
              stack.syst_only_asym->SetPointEYlow(i, std::sqrt(total_error_down * total_error_down + combined_systematics_dn.second[i] * combined_systematics_dn.second[i]));
          }
      }

      // Propagate syst errors to the stat + syst histogram
      for (uint32_t i = 1; i <= (uint32_t) stack.syst_only->GetNbinsX(); i++) {
          float syst_error = stack.syst_only->GetBinError(i);
          float stat_error = stack.stat_only->GetBinError(i);
          float syst_error_up = stack.syst_only_asym->GetErrorYhigh(i-1);
          float syst_error_dn = stack.syst_only_asym->GetErrorYlow(i-1);
          //std::cout << "syst_error " << syst_error << " stat_error " << stat_error<< " syst_error_up " << syst_error_up << " syst_error_dn " << syst_error_dn << std::endl;
          stack.stat_and_syst->SetBinError(i, std::sqrt(syst_error * syst_error + stat_error * stat_error));
          stack.stat_and_syst_asym->SetPointEYhigh(i-1, std::sqrt(syst_error_up * syst_error_up + stat_error * stat_error));
          stack.stat_and_syst_asym->SetPointEYlow(i-1, std::sqrt(syst_error_dn * syst_error_dn + stat_error * stat_error));
      }

      //////////// siglike systs
      std::map<std::string, std::vector<float>> combined_systematics_siglike_map_up;
      std::map<std::string, std::vector<float>> combined_systematics_siglike_map_dn;

      for ( auto& file: m_plotIt.getFiles([this,index] ( const File& f ) {
            return ( f.type != DATA ) && ( ! f.systematics_siglike->empty() )
                && ( ( f.type != MC ) || ( f.stack_index == index ) ) ;
            } ) ) {

          for (auto& syst: *file.systematics_siglike) {

              std::vector<float>* combined_systematics_up;
              std::vector<float>* combined_systematics_dn;
              auto map_it_up = combined_systematics_siglike_map_up.find(syst.name());
              auto map_it_dn = combined_systematics_siglike_map_dn.find(syst.name());

              if (map_it_up == combined_systematics_siglike_map_up.end()) {
                  combined_systematics_up = &combined_systematics_siglike_map_up[syst.name()];
                  combined_systematics_up->resize(stack.syst_only->GetNbinsX(), 0.);
              } else {
                  combined_systematics_up = &map_it_up->second;
              }
              if (map_it_dn == combined_systematics_siglike_map_dn.end()) {
                  combined_systematics_dn = &combined_systematics_siglike_map_dn[syst.name()];
                  combined_systematics_dn->resize(stack.syst_only->GetNbinsX(), 0.);
              } else {
                  combined_systematics_dn = &map_it_dn->second;
              }

              TH1* nominal_shape = static_cast<TH1*>(syst.nominal_shape.get());
              TH1* up_shape = static_cast<TH1*>(syst.up_shape.get());
              TH1* down_shape = static_cast<TH1*>(syst.down_shape.get());

              if (! nominal_shape || ! up_shape || ! down_shape)
                  continue;


              for (uint32_t i = 1; i <= (uint32_t) stack.syst_only->GetNbinsX(); i++) {

                  float syst_error_up = up_shape->GetBinContent(i) - nominal_shape->GetBinContent(i);
                  float syst_error_dn = down_shape->GetBinContent(i) - nominal_shape->GetBinContent(i);

                  //std::cout << "per bin " << syst.name() << " " << syst_error_up << " " << syst_error_dn << std::endl;

                  // Only propagate uncertainties for MC, not signal
                  if (file.type == MC) {
                      (*combined_systematics_up)[i - 1] += syst_error_up;
                      (*combined_systematics_dn)[i - 1] += syst_error_dn;

                      if (i == (uint32_t) stack.syst_only->GetNbinsX()) {
                          float syst_error_up_over = up_shape->GetBinContent(i+1) - nominal_shape->GetBinContent(i+1);
                          float syst_error_dn_over = down_shape->GetBinContent(i+1) - nominal_shape->GetBinContent(i+1);
                          (*combined_systematics_up)[i - 1] += syst_error_up_over;
                          (*combined_systematics_dn)[i - 1] += syst_error_dn_over;
                      }
                  }
              }
          }
      }

      // Combine all systematics in one
      // Consider that all the systematics are CORRELATED!
      for (auto& combined_systematics_up: combined_systematics_siglike_map_up) {
          for (size_t i = 1; i <= (size_t) stack.syst_only->GetNbinsX(); i++) {
              float total_error_up = stack.syst_siglike_up->GetBinContent(i);
              stack.syst_siglike_up->SetBinContent(i, total_error_up + combined_systematics_up.second[i-1]);
          }
      }

      for (auto& combined_systematics_dn: combined_systematics_siglike_map_dn) {
          for (size_t i = 1; i <= (size_t) stack.syst_only->GetNbinsX(); i++) {
              float total_error_down = stack.syst_siglike_dn->GetBinContent(i);
              stack.syst_siglike_dn->SetBinContent(i, total_error_down + combined_systematics_dn.second[i-1]);
          }
      }
      //for (size_t i = 0; i < (size_t) stack.syst_only->GetNbinsX(); i++) {
      //    std::cout << stack.syst_siglike_up->GetBinContent(i) << " " << stack.syst_siglike_dn->GetBinContent(i) << std::endl;
      //}
  }

  void TH1Plotter::computeSystematics(Stacks& stacks, Summary& summary) {
      for (auto& stack: stacks)
          computeSystematics(stack.first, stack.second, summary);
  }

  boost::optional<Summary> TH1Plotter::plot(TCanvas& c, Plot& plot) {
    c.cd();

    Summary global_summary;

    // Rescale and style histograms
    for (auto& file : m_plotIt.getFiles()) {
      setHistogramStyle(file);

      TH1* h = dynamic_cast<TH1*>(file.object);
      h->Rebin(plot.rebin);

      if (file.type != DATA) {
        plot.is_rescaled = true;

        auto generated_events = file.generated_events;

        if (m_plotIt.getConfiguration().generated_events_histogram.length() > 0 and file.generated_events < 0) {
            //generated_events = - 1.0 if not declared in file yaml
            std::shared_ptr<TFile> input(TFile::Open(file.path.c_str()));
            TH1* hevt = dynamic_cast<TH1*>(input->Get(m_plotIt.getConfiguration().generated_events_histogram.c_str()));
            generated_events = hevt->GetBinContent(m_plotIt.getConfiguration().generated_events_bin);
        }
        float factor = file.cross_section * file.branching_ratio / generated_events;

        if (! m_plotIt.getConfiguration().no_lumi_rescaling) {
          factor *= m_plotIt.getConfiguration().luminosity.at(file.era);
        }

        if (! CommandLineCfg::get().ignore_scales) {
          factor *= m_plotIt.getConfiguration().scale * file.scale;
        }

        h->Scale(factor);
        if (plot.scale_option.length() > 0)
          h->Scale(1.0, plot.scale_option.c_str());

        SummaryItem summary;
        summary.name = file.pretty_name;
        summary.process_id = file.id;

        double rescaled_integral_error = 0;
        //double rescaled_integral = h->IntegralAndError(h->GetXaxis()->GetFirst(), h->GetXaxis()->GetLast(), rescaled_integral_error);
        double rescaled_integral = h->IntegralAndError(0, h->GetNbinsX() + 1, rescaled_integral_error);

        summary.events = rescaled_integral;
        summary.events_uncertainty = rescaled_integral_error;
        summary.events_uncertainty_up = rescaled_integral_error;
        summary.events_uncertainty_down = rescaled_integral_error;

        // FIXME: Probably invalid in case of weights...

        // Bayesian efficiency
        // Taken from https://root.cern.ch/doc/master/TEfficiency_8cxx_source.html#l02428

        /*
        // Use a flat prior (equivalent to Beta(1, 1))
        float alpha = 1.;
        float beta = 1.;

        summary.efficiency = TEfficiency::BetaMean(integral + alpha, file.generated_events - integral + beta);
        summary.efficiency_uncertainty = TEfficiency::Bayesian(file.generated_events, integral, 0.68, alpha, beta, true) - summary.efficiency;
        */

        global_summary.add(file.type, summary);

        // Update all systematics for this file
        for (auto& syst: *file.systematics) {
          syst.update();

          syst.scale(factor);
          syst.rebin(plot.rebin);
          if (plot.scale_option.length() > 0)
            syst.scale(1.0, plot.scale_option.c_str());
        }
        for (auto& syst: *file.systematics_siglike) {
          syst.update();

          syst.scale(factor);
          syst.rebin(plot.rebin);
          if (plot.scale_option.length() > 0)
            syst.scale(1.0, plot.scale_option.c_str());
        }
      } else {
        SummaryItem summary;
        summary.name = file.pretty_name;
        summary.process_id = file.id;
        summary.events = h->Integral();
        global_summary.add(file.type, summary);
      }

      // Add overflow to first and last bin if requested
      if (plot.show_overflow) {
        addOverflow(h, file.type, plot);

        if (file.type != DATA) {
            for (auto& syst: *file.systematics) {
                addOverflow(static_cast<TH1*>(syst.nominal_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.up_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.down_shape.get()), file.type, plot);
            }
            for (auto& syst: *file.systematics_siglike) {
                addOverflow(static_cast<TH1*>(syst.nominal_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.up_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.down_shape.get()), file.type, plot);
            }
        }
      } else if (plot.show_onlyoverflow) {
        addOnlyOverflow(h, file.type, plot);

        if (file.type != DATA) {
            for (auto& syst: *file.systematics) {
                addOnlyOverflow(static_cast<TH1*>(syst.nominal_shape.get()), file.type, plot);
                addOnlyOverflow(static_cast<TH1*>(syst.up_shape.get()), file.type, plot);
                addOnlyOverflow(static_cast<TH1*>(syst.down_shape.get()), file.type, plot);
            }
            for (auto& syst: *file.systematics_siglike) {
                addOnlyOverflow(static_cast<TH1*>(syst.nominal_shape.get()), file.type, plot);
                addOnlyOverflow(static_cast<TH1*>(syst.up_shape.get()), file.type, plot);
                addOnlyOverflow(static_cast<TH1*>(syst.down_shape.get()), file.type, plot);
            }
        }
      }
    }

    std::shared_ptr<TH1> h_data;
    std::string data_drawing_options;

    std::vector<File> signal_files;

    for (auto& file: m_plotIt.getFiles()) {
      if (file.type == SIGNAL) {
        signal_files.push_back(file);
      } else if (file.type == DATA) {
        if (! h_data.get()) {
          h_data.reset(dynamic_cast<TH1*>(file.object->Clone()));
          h_data->SetDirectory(nullptr);
          h_data->Sumw2(false); // Disable SumW2 for data
          h_data->SetBinErrorOption((TH1::EBinErrorOpt) plot.errors_type);
          data_drawing_options += m_plotIt.getPlotStyle(file)->drawing_options;
        } else {
          h_data->Add(dynamic_cast<TH1*>(file.object));
        }
      }
    }

    auto mc_stacks = buildStacks(plot.sort_by_yields);

    if (plot.no_data || ((h_data.get()) && !h_data->GetSumOfWeights()))
      h_data.reset();

    bool has_data = h_data.get() != nullptr;
    bool has_mc = !mc_stacks.empty();

    bool no_systematics = false;

    double h_data_integral = 1.0;
    if (has_data) {
        h_data_integral = h_data->Integral();
        if (plot.scale_option.length() > 0)
            h_data->Scale(1.0, plot.scale_option.c_str());
    }

    if (plot.normalized) {
        // Normalize each plot
        for (auto& file: m_plotIt.getFiles()) {
            if (file.type == SIGNAL) {
                TH1* h = dynamic_cast<TH1*>(file.object);
                h->Scale(1. / fabs(h->GetSumOfWeights()));
            }
        }

        if (h_data.get()) {
            h_data->Scale(1. / h_data->GetSumOfWeights());
        }

        std::for_each(mc_stacks.begin(), mc_stacks.end(), [](TH1Plotter::Stacks::value_type& value) {
            TIter next(value.second.stack->GetHists());
            TH1* h = nullptr;
            while ((h = static_cast<TH1*>(next()))) {
                h->Scale(1. / std::abs(value.second.stat_only->GetSumOfWeights()));
            }

            value.second.stat_only->Scale(1. / std::abs(value.second.stat_only->GetSumOfWeights()));
        });

        no_systematics = true;
    }

    // Blind data if requested
    // It's not enough to put the bin content to zero, because
    // ROOT will show the marker, even with 'P'
    // The histogram is cloned, reset, and only the non-blinded bins are filled
    std::shared_ptr<TBox> m_blinded_area;
    if (!CommandLineCfg::get().unblind && has_data && plot.blinded_range.valid()) {
        float start = plot.blinded_range.start;
        float end = plot.blinded_range.end;

        size_t start_bin = h_data->FindBin(start);
        size_t end_bin = h_data->FindBin(end);

        TH1* clone = static_cast<TH1*>(h_data->Clone());
        clone->SetDirectory(nullptr);

        h_data->Reset();
        h_data->Sumw2(false);

        for (size_t i = 0; i < start_bin; i++) {
            h_data->SetBinContent(i, clone->GetBinContent(i));
            h_data->SetBinError(i, clone->GetBinError(i));
        }

        for (size_t i = end_bin + 1; i <= static_cast<size_t>(h_data->GetNbinsX()); i++) {
            h_data->SetBinContent(i, clone->GetBinContent(i));
            h_data->SetBinError(i, clone->GetBinError(i));
        }

        delete clone;
    }

    if (has_mc) {
        // Prepare systematics histograms
        std::for_each(mc_stacks.begin(), mc_stacks.end(), [&no_systematics](TH1Plotter::Stacks::value_type& value) {

            value.second.stat_and_syst.reset(static_cast<TH1*>(value.second.stat_only->Clone()));
            value.second.syst_siglike_up.reset(static_cast<TH1*>(value.second.stat_only->Clone()));
            value.second.syst_siglike_dn.reset(static_cast<TH1*>(value.second.stat_only->Clone()));
            value.second.stat_and_syst->SetDirectory(nullptr);
            value.second.syst_siglike_up->SetDirectory(nullptr);
            value.second.syst_siglike_dn->SetDirectory(nullptr);
            value.second.syst_siglike_up->Reset();
            value.second.syst_siglike_dn->Reset();

            uint32_t nbins = value.second.stat_only->GetNbinsX();
            float xbins[nbins];
            float xbins_err[nbins];
            float ybins[nbins];
            float zeros[nbins];

            for (uint32_t i = 0; i < nbins; i++) {
              xbins[i] = value.second.stat_only->GetBinCenter(i+1);
              xbins_err[i] = value.second.stat_only->GetBinWidth(i+1) / 2.0;
              ybins[i] = value.second.stat_only->GetBinContent(i+1);
              zeros[i] = 0.0;
            }

            auto gr_stat_and_syst = new TGraphAsymmErrors(nbins, xbins, ybins, xbins_err, xbins_err, zeros, zeros);
            value.second.stat_and_syst_asym.reset(gr_stat_and_syst);

            if (! no_systematics) {
                value.second.syst_only.reset(static_cast<TH1*>(value.second.stat_only->Clone()));
                value.second.syst_only->SetDirectory(nullptr);

                // Clear statistical errors
                for (uint32_t i = 1; i <= (uint32_t) value.second.syst_only->GetNbinsX(); i++) {
                  value.second.syst_only->SetBinError(i, 0);
                }

                auto gr_syst_only = new TGraphAsymmErrors(nbins, xbins, ybins, xbins_err, xbins_err, zeros, zeros);
                value.second.syst_only_asym.reset(gr_syst_only);
            }

        });
    }

    if (!no_systematics && plot.show_errors) {
        computeSystematics(mc_stacks, global_summary);
    }

    // Store all the histograms to draw, and find the one with the highest maximum
    std::vector<std::pair<TObject*, std::string>> toDraw = { std::make_pair(h_data.get(), data_drawing_options) };
    //for (File& signal: signal_files) {
    //  toDraw.push_back(std::make_pair(signal.object, m_plotIt.getPlotStyle(signal)->drawing_options));
    //}

    for (auto& mc_stack: mc_stacks) {
        toDraw.push_back(std::make_pair(mc_stack.second.stack.get(), ""));
    }

    // Remove NULL items
    toDraw.erase(
        std::remove_if(toDraw.begin(), toDraw.end(), [](std::pair<TObject*, std::string>& element) {
          return element.first == nullptr;
        }), toDraw.end()
      );

    if (!toDraw.size()) {
      std::cerr << "Error: nothing to draw." << std::endl;
      return boost::none;
    };

    // Sort object by minimum
    std::sort(toDraw.begin(), toDraw.end(), [&plot](std::pair<TObject*, std::string> a, std::pair<TObject*, std::string> b) {
        return (!plot.log_y) ? (getMinimum(a.first) < getMinimum(b.first)) : (getPositiveMinimum(a.first) < getPositiveMinimum(b.first));
      });

    float minimum = (!plot.log_y) ? getMinimum(toDraw[0].first) : getPositiveMinimum(toDraw[0].first);

    // Sort objects by maximum
    std::sort(toDraw.begin(), toDraw.end(), [](std::pair<TObject*, std::string> a, std::pair<TObject*, std::string> b) {
        return getMaximum(a.first) > getMaximum(b.first);
      });

    float maximum = getMaximum(toDraw[0].first);

    if (!has_data || !has_mc)
        plot.show_ratio = false;

    if (plot.show_ratio && mc_stacks.size() != 1) {
        plot.show_ratio = false;
    }

    std::shared_ptr<TPad> hi_pad;
    std::shared_ptr<TPad> low_pad;
    if (plot.show_ratio) {
      hi_pad = std::make_shared<TPad>("pad_hi", "", 0., 0.33333, 1, 1);
      hi_pad->Draw();
      const auto& config = m_plotIt.getConfiguration();
      hi_pad->SetTopMargin(config.margin_top / .6666);
      hi_pad->SetLeftMargin(config.margin_left);
      hi_pad->SetBottomMargin(0.03);
      hi_pad->SetRightMargin(config.margin_right);

      low_pad = std::make_shared<TPad>("pad_lo", "", 0., 0., 1, 0.33333);
      low_pad->Draw();
      low_pad->SetLeftMargin(config.margin_left);
      low_pad->SetTopMargin(.01);
      //low_pad->SetBottomMargin(config.margin_bottom / .3333);
      low_pad->SetBottomMargin(config.margin_bottom / .28);
      low_pad->SetRightMargin(config.margin_right);
      low_pad->SetTickx(1);

      hi_pad->cd();
      if (plot.log_y)
        hi_pad->SetLogy();

      if (plot.log_x) {
        hi_pad->SetLogx();
        low_pad->SetLogx();
      }
    }

    // Take into account systematics for maximum
    if (has_mc && !no_systematics) {
        float maximum_with_errors = 0;

        std::for_each(mc_stacks.begin(), mc_stacks.end(), [&maximum_with_errors](TH1Plotter::Stacks::value_type& value) {
            float local_max = 0;
            for (size_t b = 1; b <= (size_t) value.second.stat_and_syst->GetNbinsX(); b++) {
                float m = value.second.stat_and_syst->GetBinContent(b) + value.second.stat_and_syst->GetBinErrorUp(b);
                local_max = std::max(local_max, m);
            }

            maximum_with_errors = std::max(maximum_with_errors, local_max);
        });

        maximum = std::max(maximum, maximum_with_errors);
    }

    auto x_axis_range = plot.log_x ? plot.log_x_axis_range : plot.x_axis_range;
    auto y_axis_range = plot.log_y ? plot.log_y_axis_range : plot.y_axis_range;

    toDraw[0].first->Draw(toDraw[0].second.c_str());
    setRange(toDraw[0].first, x_axis_range, y_axis_range);

    hideTicks(toDraw[0].first, plot.x_axis_hide_ticks, plot.y_axis_hide_ticks);

    float safe_margin = .2;
    if (plot.log_y)
      safe_margin = 8;

    if (! y_axis_range.valid()) {
      maximum *= 1 + safe_margin;
      if (!plot.y_axis_auto_range) setMaximum(toDraw[0].first, maximum);
      else {
        float maxfrac = 0.45;
        float max_sig = 0.0;
        if ( signal_files.size() > 0 ) {
          std::vector<float> sigMax;
          for (File& signal: signal_files) {
            TH1* h_sig_temp = dynamic_cast<TH1*>(signal.object);
            //if (plot.signal_normalize_data and !plot.no_data) h_sig_temp->Scale(h_data->Integral()/h_sig_temp->Integral());
            if (plot.signal_normalize_data and !plot.no_data) h_sig_temp->Scale(h_data_integral/h_sig_temp->Integral());
            else if (plot.signal_normalize_data and plot.no_data) {
                auto& mc_stack_tmp = mc_stacks.begin()->second;
                h_sig_temp->Scale(mc_stack_tmp.stat_only.get()->Integral()/h_sig_temp->Integral());
            }
            sigMax.push_back(h_sig_temp->GetMaximum());
          }
          max_sig = *std::max_element(sigMax.begin(), sigMax.end());
        }
        auto& mc_stack = mc_stacks.begin()->second;
        float max_mc = 0.0;
        if (has_mc) max_mc = mc_stack.stack->GetMaximum();
        int max_data = 0.0;
        if (!plot.no_data) max_data = h_data->GetMaximum();

        if (max_mc > max_sig) {
          if (plot.log_y) {
            maxfrac = 1000;
            minimum = 0.5;
          }
          if (max_mc > max_data) setMaximum(toDraw[0].first, max_mc + max_mc*maxfrac);
          else setMaximum(toDraw[0].first, max_data + max_data*maxfrac);
        }
        else {
          if (plot.log_y) {
            maxfrac = 1000;
            minimum = 0.5;
          }
          if (max_sig > max_data) setMaximum(toDraw[0].first, max_sig + max_sig*maxfrac);
          else setMaximum(toDraw[0].first, max_data + max_data*maxfrac);
        }
      }

      if (minimum <= 0 && plot.log_y) {
        double old_minimum = minimum;
        minimum = 0.1;
        std::cout << "Warning: detected minimum is negative (" << old_minimum << ") but log scale is on. Setting minimum to " << minimum << std::endl;
      }

      if (!plot.log_y)
        minimum = minimum * (1 - std::copysign(safe_margin, minimum));

      if (plot.y_axis_show_zero && !plot.log_y)
        minimum = 0;

      setMinimum(toDraw[0].first, minimum);
    } else {
        maximum = y_axis_range.end;
        minimum = y_axis_range.start;
    }

    // First, draw MC
    if (has_mc) {

        std::for_each(mc_stacks.begin(), mc_stacks.end(), [&plot, this](TH1Plotter::Stacks::value_type& value) {
            value.second.stack->Draw("same");

            // Clear all the possible stats box remaining
            value.second.stack->GetHistogram()->SetStats(false);

            TIter next(value.second.stack->GetHists());
            TH1* h = nullptr;
            while ((h = static_cast<TH1*>(next()))) {
                h->SetStats(false);
            }

            // Then, if requested, errors
            if (plot.show_errors) {
                value.second.stat_and_syst->SetMarkerSize(0);
                value.second.stat_and_syst->SetMarkerStyle(0);
                value.second.stat_and_syst->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
                value.second.stat_and_syst->SetFillColor(m_plotIt.getConfiguration().error_fill_color);

                //value.second.stat_and_syst->Draw("E2 same");
                //TemporaryPool::get().add(value.second.stat_and_syst);

                value.second.stat_and_syst_asym->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
                value.second.stat_and_syst_asym->SetFillColor(m_plotIt.getConfiguration().error_fill_color);

                value.second.stat_and_syst_asym->Draw("2");
                TemporaryPool::get().add(value.second.stat_and_syst_asym);

                // Don't draw siglike unc in upper pane, it is not visible
                /*
                if (plot.draw_siglike_unc) {
                    std::shared_ptr<TH1> syst_siglike_up_toDraw(static_cast<TH1*>(value.second.stat_and_syst->Clone("syst_siglike_up_toDraw")));
                    std::shared_ptr<TH1> syst_siglike_dn_toDraw(static_cast<TH1*>(value.second.stat_and_syst->Clone("syst_siglike_up_toDraw")));
                    syst_siglike_up_toDraw->SetFillStyle(0);
                    syst_siglike_dn_toDraw->SetFillStyle(0);
                    syst_siglike_up_toDraw->Add(value.second.syst_siglike_up.get(), 1.0);
                    syst_siglike_dn_toDraw->Add(value.second.syst_siglike_dn.get(), 1.0);

                    syst_siglike_up_toDraw->SetLineWidth(2);
                    syst_siglike_dn_toDraw->SetLineWidth(2);
                    syst_siglike_up_toDraw->SetLineColor(3);
                    syst_siglike_dn_toDraw->SetLineColor(4);

                    syst_siglike_up_toDraw->Draw("hist same");
                    syst_siglike_dn_toDraw->Draw("hist same");
                    TemporaryPool::get().add(syst_siglike_up_toDraw);
                    TemporaryPool::get().add(syst_siglike_dn_toDraw);
                }
                */
            }
        });
    }

    // Then signal
    for (File& signal: signal_files) {
      std::string options = m_plotIt.getPlotStyle(signal)->drawing_options + " same";
      TH1* h_sig_temp = dynamic_cast<TH1*>(signal.object);
      if (plot.signal_normalize_data and !plot.no_data) {
        h_sig_temp->Scale(h_data_integral/h_sig_temp->Integral());
        if (plot.scale_option.length() > 0)
          h_sig_temp->Scale(1.0, plot.scale_option.c_str());
        h_sig_temp->Draw(options.c_str());
      }
      else if (plot.signal_normalize_data and plot.no_data) {
        auto& mc_stack_tmp = mc_stacks.begin()->second;
        h_sig_temp->Scale(mc_stack_tmp.stat_only.get()->Integral()/h_sig_temp->Integral());
        if (plot.scale_option.length() > 0)
          h_sig_temp->Scale(1.0, plot.scale_option.c_str());
        h_sig_temp->Draw(options.c_str());
      }
      else {
        if (plot.scale_option.length() > 0)
          h_sig_temp->Scale(1.0, plot.scale_option.c_str());
        h_sig_temp->Draw(options.c_str());
      }
    }

    // And finally data
    if (h_data.get()) {
      //if (plot.scale_option.length() > 0)
      //  h_data->Scale(1.0, plot.scale_option.c_str());
      data_drawing_options += " same";
      h_data->Draw(data_drawing_options.c_str());
      TemporaryPool::get().add(h_data);
    }

    // Set x and y axis titles, and default style
    for (auto& obj: toDraw) {
      setDefaultStyle(obj.first, plot, (plot.show_ratio) ? 0.6666 : 1.);
      setAxisTitles(obj.first, plot);
      hideTicks(obj.first, plot.x_axis_hide_ticks, plot.y_axis_hide_ticks);
    }

    gPad->Modified();
    gPad->Update();

    // We have the plot range. Compute the shaded area corresponding to the blinded area, if any
    if (!CommandLineCfg::get().unblind && h_data.get() && plot.blinded_range.valid()) {
        int bin_x_start = h_data->FindBin(plot.blinded_range.start);
        float x_start = h_data->GetXaxis()->GetBinLowEdge(bin_x_start);
        int bin_x_end = h_data->FindBin(plot.blinded_range.end);
        float x_end = h_data->GetXaxis()->GetBinUpEdge(bin_x_end);

        float y_start = gPad->GetUymin();
        float y_end = gPad->GetUymax();

        float lm = gPad->GetLeftMargin();
        float rm = 1. - gPad->GetRightMargin();
        //float tm = 1. - gPad->GetTopMargin();
        float bm = gPad->GetBottomMargin();

        x_start = (rm - lm) * ((x_start - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;
        x_end = (rm - lm) * ((x_end - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;

        Position legend_position = plot.legend_position;
        y_start = bm;
        y_end = legend_position.y1;

        std::string options = "NB NDC";

        if (plot.log_y) {
            //options = options + " NDC";

            //float lm = gPad->GetLeftMargin();
            //float rm = 1. - gPad->GetRightMargin();
            //float tm = 1. - gPad->GetTopMargin();
            //float bm = gPad->GetBottomMargin();

            if (plot.log_x) {
                Range x_range = getXRange(toDraw[0].first);

                x_start = (rm - lm) * ((std::log(x_start) - std::log(x_range.start)) / (std::log(x_range.end) - std::log(x_range.start))) + lm;
                x_end = (rm - lm) * ((std::log(x_end) - std::log(x_range.start)) / (std::log(x_range.end) - std::log(x_range.start))) + lm;
            }// else {
            //    x_start = (rm - lm) * ((x_start - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;
            //    x_end = (rm - lm) * ((x_end - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;
            //}

            y_start = bm;
            //y_end = tm;
            y_end = legend_position.y1;
        }

        std::shared_ptr<TPave> blinded_area(new TPave(x_start, y_start, x_end, y_end, 0, options.c_str()));
        blinded_area->SetFillStyle(m_plotIt.getConfiguration().blinded_range_fill_style);
        blinded_area->SetFillColor(m_plotIt.getConfiguration().blinded_range_fill_color);

        TemporaryPool::get().add(blinded_area);
        blinded_area->Draw("same");
    }

    auto drawLine = [&](Line& line, TVirtualPad* pad) {
        Range x_range = getXRange(toDraw[0].first);

        float y_range_start = pad->GetUymin();
        float y_range_end = pad->GetUymax();

        if (std::isnan(line.start.x))
            line.start.x = x_range.start;

        if (std::isnan(line.start.y))
            line.start.y = y_range_start;

        if (std::isnan(line.end.x))
            line.end.x = x_range.end;

        if (std::isnan(line.end.y))
            line.end.y = y_range_end;

        std::shared_ptr<TLine> l(new TLine(line.start.x, line.start.y, line.end.x, line.end.y));
        TemporaryPool::get().add(l);

        l->SetLineColor(line.style->line_color);
        l->SetLineWidth(line.style->line_width);
        l->SetLineStyle(line.style->line_type);

        l->Draw("same");
    };

    for (Line& line: plot.lines) {
      // Only keep TOP lines
      if (line.pad != TOP)
        continue;

      drawLine(line, hi_pad ? hi_pad.get() : gPad);
    }

    // Redraw only axis
    toDraw[0].first->Draw("axis same");

    if (plot.show_ratio) {

      // Compute ratio and draw it
      low_pad->cd();
      low_pad->SetGridy();

      std::shared_ptr<TH1> h_low_pad_axis(static_cast<TH1*>(h_data->Clone()));
      h_low_pad_axis->SetDirectory(nullptr);
      h_low_pad_axis->Reset(); // Keep binning
      setRange(h_low_pad_axis.get(), x_axis_range, plot.ratio_y_axis_range);

      if (gROOT->GetVersionInt() >= 62204)
        setDefaultStyle(h_low_pad_axis.get(), plot, 1.);
      else
        setDefaultStyle(h_low_pad_axis.get(), plot, 3.); //old ROOT
      h_low_pad_axis->GetYaxis()->SetTitle(plot.ratio_y_axis_title.c_str());
      h_low_pad_axis->GetYaxis()->SetTickLength(0.04);
      h_low_pad_axis->GetYaxis()->SetNdivisions(505, true);
      h_low_pad_axis->GetXaxis()->SetTickLength(0.07);

      hideTicks(h_low_pad_axis.get(), plot.x_axis_hide_ticks, plot.y_axis_hide_ticks);

      h_low_pad_axis->Draw();

      auto& mc_stack = mc_stacks.begin()->second;

      //std::shared_ptr<TGraphAsymmErrors> ratio = getRatio(h_data.get(), mc_stack.stat_only.get());
      std::shared_ptr<TGraphAsymmErrors> ratio = getRatio2(h_data.get(), mc_stack.stat_only.get());
      ratio->Draw((m_plotIt.getConfiguration().ratio_style + "same").c_str());

      if (plot.ratio_y_axis_auto_range) {
        float ratio_max = plot.ratio_y_axis_range.end;
        float ratio_min = plot.ratio_y_axis_range.start;
        //float current_max = TMath::MaxElement(ratio->GetN(), ratio->GetY());
        //float current_min = TMath::MinElement(ratio->GetN(), ratio->GetY());
        auto hratio = ratio->GetHistogram();
        float current_max = hratio->GetMaximum() + hratio->GetBinError(hratio->GetMaximumBin());
        float current_min = hratio->GetMinimum() - hratio->GetBinError(hratio->GetMinimumBin());

        if(current_max > ratio_max)
          ratio_max = current_max * 1.001;
        if(current_min < ratio_min)
          ratio_min = current_min * 0.999;

        Range ratio_auto_range = {ratio_min, ratio_max};
        setRange(h_low_pad_axis.get(), x_axis_range, ratio_auto_range);
        h_low_pad_axis->Draw();
      }

      // Compute stat errors
      std::shared_ptr<TH1> h_mcstat(static_cast<TH1*>(h_low_pad_axis->Clone()));
      h_mcstat->SetDirectory(nullptr);
      h_mcstat->Reset(); // Keep binning
      h_mcstat->SetMarkerSize(0);
      h_mcstat->SetBinErrorOption(TH1::kNormal);

      if (plot.ratio_draw_mcstat_error) {

        for (uint32_t i = 1; i <= (uint32_t) h_mcstat->GetNbinsX(); i++) {

          if (mc_stack.stat_only->GetBinContent(i) == 0 || mc_stack.stat_only->GetBinError(i) == 0)
            continue;

          // relative error, delta X / X
          float stat = mc_stack.stat_only->GetBinError(i) / mc_stack.stat_only->GetBinContent(i);

          h_mcstat->SetBinContent(i, 1);
          h_mcstat->SetBinError(i, stat);
        }

        h_mcstat->SetFillStyle(m_plotIt.getConfiguration().staterror_fill_style);
        h_mcstat->SetFillColor(m_plotIt.getConfiguration().staterror_fill_color);
        setRange(h_mcstat.get(), x_axis_range, {});
        h_mcstat->Draw("E2 same");

        auto mc_legend = new TLegend(0.20, 0.82, 0.35, 0.92);
        mc_legend->AddEntry(h_mcstat.get(), "Stat. Unc.", "F");
        mc_legend->SetTextFont(42);
        mc_legend->SetFillStyle(0);
        mc_legend->SetBorderSize(0);
        mc_legend->Draw();
      }

      h_low_pad_axis->Draw("same");


      // Compute systematic errors
      std::shared_ptr<TH1> h_systematics(static_cast<TH1*>(h_low_pad_axis->Clone()));
      h_systematics->SetDirectory(nullptr);
      h_systematics->Reset(); // Keep binning
      h_systematics->SetMarkerSize(0);
      // Workaround a bug introduced in ROOT 6.08 causing the bin error to not reset to normal
      // See https://sft.its.cern.ch/jira/browse/ROOT-8808 for more details
      h_systematics->SetBinErrorOption(TH1::kNormal);

      std::shared_ptr<TH1> h_syst_siglike_up(static_cast<TH1*>(h_systematics->Clone()));
      std::shared_ptr<TH1> h_syst_siglike_dn(static_cast<TH1*>(h_systematics->Clone()));
      h_syst_siglike_up->SetDirectory(nullptr);
      h_syst_siglike_dn->SetDirectory(nullptr);
      h_syst_siglike_up->Reset();
      h_syst_siglike_dn->Reset();

      // Asymmetric unc band
      uint32_t nbins = h_systematics->GetNbinsX();
      float xbins[nbins];
      float xbins_err[nbins];
      float ybins[nbins];
      float yerrup[nbins];
      float yerrdn[nbins];

      for (uint32_t i = 0; i < nbins; i++) {
        xbins[i] = h_data->GetBinCenter(i+1);
        xbins_err[i] = h_data->GetBinWidth(i+1) / 2.0;
        ybins[i] = 1.0;
        yerrup[i] = 0.0;
        yerrdn[i] = 0.0;
      }

      bool has_syst = false;
      if (! no_systematics) {
        for (uint32_t i = 1; i <= (uint32_t) h_systematics->GetNbinsX(); i++) {

          if (mc_stack.syst_only->GetBinContent(i) == 0)
            continue;

          const auto& config = m_plotIt.getConfiguration();

          if (mc_stack.syst_only->GetBinError(i) != 0) {
            // relative error, delta X / X
            // Post-fit unc is computed with stat+syst together from the fit.
            // MC stat unc must be removed, total postfit unc includes it.
            float syst = 0.;
            if (config.syst_only or plot.post_fit) syst = mc_stack.syst_only->GetBinError(i) / mc_stack.syst_only->GetBinContent(i);
            else syst = mc_stack.stat_and_syst->GetBinError(i) / mc_stack.syst_only->GetBinContent(i);

            //std::cout << mc_stack.stat_and_syst->GetBinError(i) - mc_stack.syst_only->GetBinError(i) << std::endl;
            //std::cout << "old    " << mc_stack.stat_and_syst->GetBinError(i) << "/" << mc_stack.syst_only->GetBinContent(i) << std::endl;
            //std::cout << "new up " << mc_stack.stat_and_syst_asym->GetErrorYhigh(i) << "/" << mc_stack.syst_only->GetBinContent(i) << std::endl;
            //std::cout << "new dn " << mc_stack.stat_and_syst_asym->GetErrorYlow(i) << "/" << mc_stack.syst_only->GetBinContent(i) << std::endl;

            h_systematics->SetBinContent(i, 1);
            h_systematics->SetBinError(i, syst);
            //if (h_data.get()->GetBinContent(i) > 0) h_systematics->SetBinError(i, syst);

            has_syst = true;
          }

          if (mc_stack.syst_only_asym->GetErrorYhigh(i-1) != 0 or mc_stack.syst_only_asym->GetErrorYlow(i-1) != 0) {
            if (config.syst_only or plot.post_fit) {
              yerrup[i-1] = mc_stack.syst_only_asym->GetErrorYhigh(i-1) / mc_stack.syst_only->GetBinContent(i);
              yerrdn[i-1] = mc_stack.syst_only_asym->GetErrorYlow(i-1) / mc_stack.syst_only->GetBinContent(i);
            } else {
              yerrup[i-1] = mc_stack.stat_and_syst_asym->GetErrorYhigh(i-1) / mc_stack.syst_only->GetBinContent(i);
              yerrdn[i-1] = mc_stack.stat_and_syst_asym->GetErrorYlow(i-1) / mc_stack.syst_only->GetBinContent(i);
            }
            has_syst = true;
          }

          for (uint32_t i = 1; i <= (uint32_t) h_systematics->GetNbinsX(); i++) {

              if (mc_stack.syst_siglike_up->GetBinContent(i) != 0) {
                  float syst = (mc_stack.syst_siglike_up->GetBinContent(i) + mc_stack.syst_only->GetBinContent(i)) / mc_stack.syst_only->GetBinContent(i);
                  h_syst_siglike_up->SetBinContent(i, syst);
                  h_syst_siglike_up->SetBinError(i, 0);
              }
              if (mc_stack.syst_siglike_dn->GetBinContent(i) != 0) {
                  float syst = (mc_stack.syst_siglike_dn->GetBinContent(i) + mc_stack.syst_only->GetBinContent(i)) / mc_stack.syst_only->GetBinContent(i);
                  h_syst_siglike_dn->SetBinContent(i, syst);
                  h_syst_siglike_dn->SetBinError(i, 0);
              }
          }
        }
      }

      std::shared_ptr<TGraphAsymmErrors> graph_systematics(new TGraphAsymmErrors(nbins, xbins, ybins, xbins_err, xbins_err, yerrdn, yerrup));

      if (has_syst) {
        h_systematics->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
        h_systematics->SetFillColor(m_plotIt.getConfiguration().error_fill_color);
        setRange(h_systematics.get(), x_axis_range, {});
        //h_systematics->Draw("E2 same");

        graph_systematics->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
        graph_systematics->SetFillColor(m_plotIt.getConfiguration().error_fill_color);
        setRange(graph_systematics.get(), x_axis_range, {});
        graph_systematics->Draw("2");


        if (plot.draw_siglike_unc) {
            h_syst_siglike_up->SetFillStyle(0);
            h_syst_siglike_dn->SetFillStyle(0);
            h_syst_siglike_up->SetLineWidth(3);
            h_syst_siglike_dn->SetLineWidth(3);
            h_syst_siglike_up->SetLineColor(601);
            h_syst_siglike_dn->SetLineColor(633);
            h_syst_siglike_up->Draw("same hist");
            h_syst_siglike_dn->Draw("same hist");

            auto syst_legend = new TLegend(0.75, 0.45, 0.95, 0.55);
            syst_legend->SetNColumns(2);
            syst_legend->AddEntry(h_syst_siglike_up.get(), "Var. up", "F");
            syst_legend->AddEntry(h_syst_siglike_dn.get(), "Var. dn", "F");
            syst_legend->SetTextFont(42);
            syst_legend->SetFillStyle(0);
            syst_legend->SetBorderSize(0);
            syst_legend->Draw();

        }
      }

      h_low_pad_axis->Draw("same");

      if (plot.fit_ratio) {
        float xMin, xMax;
        if (plot.ratio_fit_range.valid()) {
          xMin = plot.ratio_fit_range.start;
          xMax = plot.ratio_fit_range.end;
        } else {
          xMin = h_low_pad_axis->GetXaxis()->GetBinLowEdge(1);
          xMax = h_low_pad_axis->GetXaxis()->GetBinUpEdge(h_low_pad_axis->GetXaxis()->GetLast());
        }

        std::shared_ptr<TF1> fct = std::make_shared<TF1>("fit_function", plot.ratio_fit_function.c_str(), xMin, xMax);
        fct->SetNpx(m_plotIt.getConfiguration().ratio_fit_n_points);

        TFitResultPtr fit_result = ratio->Fit(fct.get(), "SMRNEQ");
        if (fit_result->IsValid()) {
          std::shared_ptr<TH1> errors = std::make_shared<TH1D>("errors", "errors", m_plotIt.getConfiguration().ratio_fit_n_points, xMin, xMax);
          errors->SetDirectory(nullptr);
          (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors.get(), 0.68);
          errors->SetStats(false);
          errors->SetMarkerSize(0);
          errors->SetFillColor(m_plotIt.getConfiguration().ratio_fit_error_fill_color);
          errors->SetFillStyle(m_plotIt.getConfiguration().ratio_fit_error_fill_style);
          errors->Draw("e3 same");

          fct->SetLineWidth(m_plotIt.getConfiguration().ratio_fit_line_width);
          fct->SetLineColor(m_plotIt.getConfiguration().ratio_fit_line_color);
          fct->SetLineStyle(m_plotIt.getConfiguration().ratio_fit_line_style);
          fct->Draw("same");

          if (plot.ratio_fit_legend.length() > 0) {
            uint32_t fit_parameters = fct->GetNpar();
            boost::format formatter = get_formatter(plot.ratio_fit_legend);

            for (uint32_t i = 0; i < fit_parameters; i++) {
              formatter % fct->GetParameter(i);
            }

            std::string legend = formatter.str();

            std::shared_ptr<TLatex> t(new TLatex(plot.ratio_fit_legend_position.x, plot.ratio_fit_legend_position.y, legend.c_str()));
            t->SetNDC(true);
            t->SetTextFont(62);
            //t->SetTextSize(LABEL_FONTSIZE - 4);
            t->SetTextSize(0.65*0.2);
            t->Draw();

            TemporaryPool::get().add(t);
          }

          TemporaryPool::get().add(errors);
          TemporaryPool::get().add(fct);
        }
      }

      h_low_pad_axis->Draw("same");
      ratio->Draw((m_plotIt.getConfiguration().ratio_style + "same").c_str());

      // Hide top pad label
      hideXTitle(toDraw[0].first);

      low_pad->Modified();
      low_pad->Update();

      for (Line& line: plot.lines) {
        // Only keep BOTTOM lines
        if (line.pad != BOTTOM)
          continue;

        drawLine(line, low_pad.get());
      }

      TemporaryPool::get().add(h_low_pad_axis);
      TemporaryPool::get().add(ratio);
      //TemporaryPool::get().add(h_systematics);
      TemporaryPool::get().add(h_mcstat);
      TemporaryPool::get().add(graph_systematics);
      TemporaryPool::get().add(hi_pad);
      TemporaryPool::get().add(low_pad);
      if (plot.draw_siglike_unc) {
          TemporaryPool::get().add(h_syst_siglike_up);
          TemporaryPool::get().add(h_syst_siglike_dn);
      }
    }

    if (has_mc && mc_stacks.size() == 1 && plot.fit) {

      auto& mc_stack = mc_stacks.begin()->second;

      float xMin, xMax;
      if (plot.fit_range.valid()) {
        xMin = plot.fit_range.start;
        xMax = plot.fit_range.end;
      } else {
        xMin = mc_stack.stat_only->GetXaxis()->GetBinLowEdge(1);
        xMax = mc_stack.stat_only->GetXaxis()->GetBinUpEdge(mc_stack.stat_only->GetXaxis()->GetLast());
      }

      std::shared_ptr<TF1> fct = std::make_shared<TF1>("fit_function", plot.fit_function.c_str(), xMin, xMax);
      fct->SetNpx(m_plotIt.getConfiguration().fit_n_points);

      TH1* mc_hist = mc_stack.stat_only.get();
      TFitResultPtr fit_result = mc_hist->Fit(fct.get(), "SMRNEQ");
      if (fit_result->IsValid()) {
        std::shared_ptr<TH1> errors = std::make_shared<TH1D>("errors", "errors", m_plotIt.getConfiguration().fit_n_points, xMin, xMax);
        errors->SetDirectory(nullptr);
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors.get(), 0.68);
        errors->SetStats(false);
        errors->SetMarkerSize(0);
        errors->SetFillColor(m_plotIt.getConfiguration().fit_error_fill_color);
        errors->SetFillStyle(m_plotIt.getConfiguration().fit_error_fill_style);
        errors->Draw("e3 same");

        fct->SetLineWidth(m_plotIt.getConfiguration().fit_line_width);
        fct->SetLineColor(m_plotIt.getConfiguration().fit_line_color);
        fct->SetLineStyle(m_plotIt.getConfiguration().fit_line_style);
        fct->Draw("same");

        if (plot.fit_legend.length() > 0) {
          uint32_t fit_parameters = fct->GetNpar();
          boost::format formatter = get_formatter(plot.fit_legend);

          for (uint32_t i = 0; i < fit_parameters; i++) {
            formatter % fct->GetParameter(i);
          }

          std::string legend = formatter.str();

          std::shared_ptr<TLatex> t(new TLatex(plot.fit_legend_position.x, plot.fit_legend_position.y, legend.c_str()));
          t->SetNDC(true);
          t->SetTextFont(62);
          t->SetTextSize(LABEL_FONTSIZE - 0);
          t->Draw();

          TemporaryPool::get().add(t);
        }

        TemporaryPool::get().add(errors);
        TemporaryPool::get().add(fct);
      }
    }

    gPad->Modified();
    gPad->Update();
    gPad->RedrawAxis();

    if (hi_pad.get())
      hi_pad->cd();

    return global_summary;
  }

  void TH1Plotter::setHistogramStyle(const File& file) {
    TH1* h = dynamic_cast<TH1*>(file.object);

    std::shared_ptr<PlotStyle> style = m_plotIt.getPlotStyle(file);

    if (style->fill_color != -1)
      h->SetFillColor(style->fill_color);

    if (style->fill_type != -1)
      h->SetFillStyle(style->fill_type);

    if (style->line_color != -1)
      h->SetLineColor(style->line_color);

    if (style->line_width != -1)
      h->SetLineWidth(style->line_width);

    if (style->line_type != -1)
      h->SetLineStyle(style->line_type);

    if (style->marker_size != -1)
      h->SetMarkerSize(style->marker_size);

    if (style->marker_color != -1)
      h->SetMarkerColor(style->marker_color);

    if (style->marker_type != -1)
      h->SetMarkerStyle(style->marker_type);

    if (file.type == MC && style->line_color == -1 && style->fill_color != -1)
      h->SetLineColor(style->fill_color);
  }

  void TH1Plotter::addOverflow(TH1* h, Type type, const Plot& plot) {

    if (!h || !h->GetEntries())
        return;

    size_t first_bin = 1;
    size_t last_bin = h->GetNbinsX();

    auto x_axis_range = plot.log_x ? plot.log_x_axis_range : plot.x_axis_range;

    if (x_axis_range.valid()) {
      std::shared_ptr<TH1> copy(dynamic_cast<TH1*>(h->Clone()));
      copy->SetDirectory(nullptr);
      copy->GetXaxis()->SetRangeUser(x_axis_range.start, x_axis_range.end);

      // Find first and last bin corresponding to the given range
      first_bin = copy->GetXaxis()->GetFirst();
      last_bin = copy->GetXaxis()->GetLast();
    }

    // GetBinError returns sqrt(SumW2) for a given bin
    // SetBinError updates SumW2 for a given bin with error*error

    float underflow = 0;
    float underflow_sumw2 = 0;
    for (size_t i = 0; i < first_bin; i++) {
      underflow += h->GetBinContent(i);
      underflow_sumw2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    float overflow = 0;
    float overflow_sumw2 = 0;
    for (size_t i = last_bin + 1; i <= (size_t) h->GetNbinsX() + 1; i++) {
      overflow += h->GetBinContent(i);
      overflow_sumw2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    // Clear out-of-range bin content so that Integral() still returns the right value
    for (size_t i = 1; i < first_bin; i++) {
      h->SetBinContent(i, 0);
    }
    for (size_t i = last_bin + 1; i < (size_t) h->GetNbinsX() + 1; i++) {
      h->SetBinContent(i, 0);
    }
    // Clear also underflow and overflow bins (SetBinContent on these may try to extend the axes)
    h->ClearUnderflowAndOverflow();

    float first_bin_content = h->GetBinContent(first_bin);
    float first_bin_sumw2 = h->GetBinError(first_bin) * h->GetBinError(first_bin);

    float last_bin_content = h->GetBinContent(last_bin);
    float last_bin_sumw2 = h->GetBinError(last_bin) * h->GetBinError(last_bin);

    h->SetBinContent(first_bin, first_bin_content + underflow);
    if (type != DATA)
        h->SetBinError(first_bin, sqrt(underflow_sumw2 + first_bin_sumw2));

    h->SetBinContent(last_bin, last_bin_content + overflow);
    if (type != DATA)
        h->SetBinError(last_bin, sqrt(overflow_sumw2 + last_bin_sumw2));
  }

  void TH1Plotter::addOnlyOverflow(TH1* h, Type type, const Plot& plot) {

    if (!h || !h->GetEntries())
        return;

    size_t last_bin = h->GetNbinsX();

    auto x_axis_range = plot.log_x ? plot.log_x_axis_range : plot.x_axis_range;

    if (x_axis_range.valid()) {
      std::shared_ptr<TH1> copy(dynamic_cast<TH1*>(h->Clone()));
      copy->SetDirectory(nullptr);
      copy->GetXaxis()->SetRangeUser(x_axis_range.start, x_axis_range.end);

      // Find first and last bin corresponding to the given range
      last_bin = copy->GetXaxis()->GetLast();
    }

    // GetBinError returns sqrt(SumW2) for a given bin
    // SetBinError updates SumW2 for a given bin with error*error

    float overflow = 0;
    float overflow_sumw2 = 0;
    for (size_t i = last_bin + 1; i <= (size_t) h->GetNbinsX() + 1; i++) {
      overflow += h->GetBinContent(i);
      overflow_sumw2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    // Clear out-of-range bin content so that Integral() still returns the right value
    for (size_t i = last_bin + 1; i < (size_t) h->GetNbinsX() + 1; i++) {
      h->SetBinContent(i, 0);
    }
    // Clear also underflow and overflow bins (SetBinContent on these may try to extend the axes)
    h->ClearUnderflowAndOverflow();

    float last_bin_content = h->GetBinContent(last_bin);
    float last_bin_sumw2 = h->GetBinError(last_bin) * h->GetBinError(last_bin);

    h->SetBinContent(last_bin, last_bin_content + overflow);
    if (type != DATA)
        h->SetBinError(last_bin, sqrt(overflow_sumw2 + last_bin_sumw2));
  }
}
