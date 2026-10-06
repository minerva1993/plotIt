#pragma once

#include <plotter.h>

namespace plotIt {
    class TH1Plotter: public plotter {
        public:
            struct Stack {
                std::shared_ptr<THStack> stack;
                std::shared_ptr<TH1> stat_only;
                std::shared_ptr<TH1> syst_only;
                std::shared_ptr<TH1> stat_and_syst;
                std::shared_ptr<TGraphAsymmErrors> syst_only_asym;
                std::shared_ptr<TGraphAsymmErrors> stat_and_syst_asym;
                std::shared_ptr<TH1> syst_siglike_up;
                std::shared_ptr<TH1> syst_siglike_dn;
                std::shared_ptr<TH1> syst_mass_up;
                std::shared_ptr<TH1> syst_mass_dn;
                std::shared_ptr<TH1> syst_width_up;
                std::shared_ptr<TH1> syst_width_dn;
            };

            using Stacks = std::vector<std::pair<int64_t, Stack>>;

            TH1Plotter(plotIt& plotIt):
                plotter(plotIt) {
                }

            virtual boost::optional<Summary> plot(TCanvas& c, Plot& plot);
            virtual bool supports(TObject& object);

        private:
            void setHistogramStyle(const File& file);
            void addOverflow(TH1* h, Type type, const Plot& plot);
            void addOnlyOverflow(TH1* h, Type type, const Plot& plot);

            Stack buildStack(int64_t index, bool sortByYields);
            Stacks buildStacks(bool sortByYields);

            void computeSystematics(int64_t index, Stack& stack, Summary& summary);
            void computeSystematics(Stacks& stacks, Summary& summary);
    };
}
