#pragma once

#include <boost/filesystem.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <memory>
#include <iomanip>
#include <iostream>

#include "yaml-cpp/yaml.h"

#include <TH1.h>
#include <THStack.h>
#include <TStyle.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>

#include <vector>
#include <string>
#include <glob.h>
#include <unordered_map>

#include <types.h>
#include <defines.h>
#include <uuid.h>

namespace YAML {
  class Node;
}

class TFile;
class TObject;
class TCanvas;
class TLegend;

namespace fs = boost::filesystem;

namespace plotIt {
  
  class plotIt {
    public:
      plotIt(const fs::path& outputPath);
      bool parseConfigurationFile(const std::string& file, const fs::path& histogramsPath);
      void plotAll();

      // a bit of infrastructure to retrieve selected file lists
      // stored as vector<const*> but behaving as reference vectors
      class file_list {
      public:
        file_list() = default;
        std::size_t size() const { return m_sFiles.size(); }
        void push_back(const File& file) { m_sFiles.push_back(&file); }
        const File& operator[] ( std::size_t idx ) const { return *(m_sFiles[idx]); }
        class const_iterator : public boost::iterator_adaptor<const_iterator,
            std::vector<const File*>::const_iterator,File,
            boost::random_access_traversal_tag,const File&> {
          public:
            explicit const_iterator(const base_type& base) : iterator_adaptor_(base) {}
            reference dereference() const { return **base_reference(); }
        };
        using value_type = const_iterator::value_type;
        const_iterator begin() const { return const_iterator{std::begin(m_sFiles)}; }
        const_iterator end() const { return const_iterator{std::end(m_sFiles)}; }
      private:
        std::vector<const File*> m_sFiles;
      };
      template<typename Predicate>
      file_list getFiles(Predicate pred) const {
        file_list result;
        std::copy_if(std::begin(m_files), std::end(m_files), std::back_inserter(result), pred);
        return result;
      }

      bool filter_eras(const File& file) const {
        return m_config.eras.empty() || file.era.empty() || ( std::end(m_config.eras) != std::find(std::begin(m_config.eras), std::end(m_config.eras), file.era) );
      }
      file_list getFiles() const { return getFiles(std::bind(&plotIt::filter_eras, this, std::placeholders::_1)); }

      const Configuration& getConfiguration() const {
        return m_config;
      }

      std::shared_ptr<PlotStyle> getPlotStyle(const File& file);

      friend PlotStyle;

    private:
      void checkOrThrow(YAML::Node& node, const std::string& name, const std::string& file);
      void parseIncludes(YAML::Node& node, const fs::path& base);
      void parseSystematicsNode(const YAML::Node& node);
      void parseFileNode(File& file, const YAML::Node& key, const YAML::Node& value);
      void parseFileNode(File& file, const YAML::Node& node);

      // Plot method
      bool plot(Plot& plot);
      bool yields(std::vector<Plot>::iterator plots_begin, std::vector<Plot>::iterator plots_end);
      bool systematics(std::vector<Plot>::iterator plots_begin, std::vector<Plot>::iterator plots_end);

      bool expandFiles();
      bool expandObjects(File& file, std::vector<Plot>& plots);
      bool loadAllObjects(File& file, std::vector<Plot>::const_iterator plots_begin, std::vector<Plot>::const_iterator plots_end);
      bool loadObject(File& file, const Plot& plot);

      void fillLegend(TLegend& legend, const Plot& plot, bool with_uncertainties);

      void parseLumiLabel();

      std::vector<Label> mergeLabels(const std::vector<Label>& labels);

      fs::path m_outputPath;

      std::vector<File> m_files;
      std::vector<Plot> m_plots;
      std::vector<SystematicPtr> m_systematics;
      std::vector<SystematicPtr> m_systematics_siglike;
      std::vector<SystematicPtr> m_systematics_mass;
      std::vector<SystematicPtr> m_systematics_width;
      std::map<std::string, Group> m_legend_groups;
      std::map<std::string, Group> m_yields_groups;

      std::unordered_map<std::string, TDirectory*> m_book_keeping_folders;

      // Current style
      std::shared_ptr<TStyle> m_style;

      Legend m_legend;
      Configuration m_config;
  };
};
