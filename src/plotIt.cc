#include "plotIt.h"

// For fnmatch()
#include <fnmatch.h>

#include <TROOT.h>
#include <TList.h>
#include <TCollection.h>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TGaxis.h>
#include <Math/QuantFuncMathCore.h>

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>

#include "tclap/CmdLine.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <commandlinecfg.h>
#include <plotters.h>
#include <pool.h>
#include <summary.h>
#include <systematics.h>
#include <utilities.h>


namespace fs = boost::filesystem;
using std::setw;

// Load libdpm at startup, on order to be sure that rfio files are working
#include <dlfcn.h>
struct Dummy
{
  Dummy()
  {
    dlopen("libdpm.so", RTLD_NOW|RTLD_GLOBAL);
  }
};
static Dummy foo;

namespace plotIt {

  plotIt::plotIt(const fs::path& outputPath):
    m_outputPath(outputPath) {

      createPlotters(*this);

      gErrorIgnoreLevel = kError;

      TH1::AddDirectory(false);
    }

  // Replace the "include" fields by the content they point to
  void plotIt::parseIncludes(YAML::Node& node, const fs::path& base) {


    if (node["include"]) {
        std::vector<std::string> files = node["include"].as<std::vector<std::string>>();

        YAML::Node merged_node;

        for (std::string& file: files) {
          fs::path ifp = fs::absolute(fs::path(file), base);
          YAML::Node root;
          try {
            root = YAML::LoadFile(ifp.string());
          } catch ( const YAML::BadFile& e ) {
            std::cout << "Problem parsing YAML file '" << ifp << "'" << std::endl;
            throw e;
          }
          if (root.Type() == YAML::NodeType::Map) {
            parseIncludes(root, ifp.parent_path());
          }

          for (YAML::const_iterator it = root.begin(); it != root.end(); ++it) {
            if (root.Type() == YAML::NodeType::Map) {
                merged_node[it->first.as<std::string>()] = it->second;
            } else {
                merged_node.push_back(*it);
            }
          }
        }

        node = merged_node;

    }

    for (YAML::iterator it = node.begin(); it != node.end(); ++it) {
      if (it->second.Type() == YAML::NodeType::Map) {
          parseIncludes(it->second, base);
      }
    }

  }

  void plotIt::parseSystematicsNode(const YAML::Node& node) {

      std::string type;
      std::string name;
      YAML::Node configuration;

      switch (node.Type()) {
        case YAML::NodeType::Scalar:
            name = node.as<std::string>();
            type = "shape";
            break;

        case YAML::NodeType::Map: {
            const auto& it = *node.begin();

            if (it.second.IsScalar())
                type = "const";
            else if (it.second["type"])
                type = it.second["type"].as<std::string>();
            else
                type = "shape";

            name = it.first.as<std::string>();
            configuration = it.second;
            } break;

        default:
            throw YAML::ParserException(node.Mark(), "Invalid systematics node. Must be either a string or a map");
      }

      if (type == "siglike") m_systematics_siglike.push_back(SystematicFactory::create(name, type, configuration));
      else m_systematics.push_back(SystematicFactory::create(name, type, configuration));
  }

  std::vector<RenameOp> parseRenameNode(const YAML::Node& node) {
      std::vector<RenameOp> ops;

      if (! node["rename"])
          return ops;

      const auto& rename_node = node["rename"];

      for (YAML::const_iterator it = rename_node.begin(); it != rename_node.end(); ++it) {
          const YAML::Node& rename_op_node = *it;
          RenameOp op;
          op.from = std::regex(rename_op_node["from"].as<std::string>(), std::regex::extended);
          op.to = rename_op_node["to"].as<std::string>();

          ops.push_back(op);
      }

      return ops;
  }

  void plotIt::parseFileNode(File& file, const YAML::Node& key, const YAML::Node& value) {

      file.path = key.as<std::string>();
      parseFileNode(file, value);
  }

  void plotIt::parseFileNode(File& file, const YAML::Node& node) {

      if (node["file"]) {
          file.path = node["file"].as<std::string>();
      }

      // Normalize path
      fs::path root = fs::path(m_config.root);
      fs::path path = fs::path(file.path);
      file.path = (root / path).string();

      if (node["pretty-name"]) {
        file.pretty_name = node["pretty-name"].as<std::string>();
      } else {
        file.pretty_name = path.stem().native();
      }

      if (node["era"])
        file.era = node["era"].as<std::string>();

      if (node["type"]) {
        std::string type = node["type"].as<std::string>();
        file.type = string_to_type(type);
      }

      if (node["scale"])
        file.scale = node["scale"].as<float>();

      if (node["cross-section"])
        file.cross_section = node["cross-section"].as<float>();

      if (node["branching-ratio"])
        file.branching_ratio = node["branching-ratio"].as<float>();

      if (node["generated-events"])
        file.generated_events = node["generated-events"].as<float>();

      if (node["order"])
        file.order = node["order"].as<int16_t>();

      if (node["group"])
        file.legend_group = node["group"].as<std::string>();

      if (node["yields-group"])
        file.yields_group = node["yields-group"].as<std::string>();

      if (node["stack-index"])
        file.stack_index = node["stack-index"].as<int64_t>();

      file.renaming_ops = parseRenameNode(node);

      file.plot_style = std::make_shared<PlotStyle>();
      file.plot_style->loadFromYAML(node, file.type);
  }

  bool plotIt::parseConfigurationFile(const std::string& file, const fs::path& histogramsPath) {
    YAML::Node f;
    try {
      f = YAML::LoadFile(file);
    } catch ( const YAML::BadFile& e ) {
      std::cout << "Problem parsing YAML file '" << file << "'" << std::endl;
      throw e;
    }

    if (CommandLineCfg::get().verbose) {
        std::cout << "Parsing configuration file ...";
    }

    parseIncludes(f, fs::absolute(fs::path(file)).parent_path());

    if (! f["files"]) {
      throw YAML::ParserException(YAML::Mark::null_mark(), "Your configuration file must have a 'files' list");
    }

    const auto& parseLabelsNode = [](YAML::Node& node) -> std::vector<Label> {
      std::vector<Label> labels;

      for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
        const YAML::Node& labelNode = *it;

        Label label;
        label.text = labelNode["text"].as<std::string>();
        label.position = labelNode["position"].as<Point>();

        if (labelNode["size"])
//          label.size = labelNode["size"].as<uint32_t>();
          label.size = labelNode["size"].as<float>();

        if (labelNode["font"])
          label.font = labelNode["font"].as<int>();

        labels.push_back(label);
      }

      return labels;
    };

    // Retrieve legend configuration
    if (f["legend"]) {
      YAML::Node node = f["legend"];

      if (node["position"])
        m_legend.position = node["position"].as<Position>();

      if (node["columns"])
          m_legend.columns = node["columns"].as<size_t>();

      if (node["entries"]) {
        const YAML::Node& node_entries = node["entries"];
        for (const auto& node_entry: node_entries) {
            std::string label = node_entry["label"].as<std::string>();

            Type type = MC;
            if (node_entry["type"])
                type = string_to_type(node_entry["type"].as<std::string>());

            int16_t order = 0;
            if (node_entry["order"])
                order = node_entry["order"].as<int16_t>();

            m_config.static_legend_entries[type].push_back(LegendEntry(nullptr, label, "", order));
        }
      }
    }


    // Retrieve plotIt configuration
    if (f["configuration"]) {
      YAML::Node node = f["configuration"];

      if (node["width"])
        m_config.width = node["width"].as<float>();

      if (node["height"])
        m_config.height = node["height"].as<float>();

      if (node["margin-left"])
        m_config.margin_left = node["margin-left"].as<float>();

      if (node["margin-right"])
        m_config.margin_right = node["margin-right"].as<float>();

      if (node["margin-top"])
        m_config.margin_top = node["margin-top"].as<float>();

      if (node["margin-bottom"])
        m_config.margin_bottom = node["margin-bottom"].as<float>();

      if (node["experiment"])
        m_config.experiment = node["experiment"].as<std::string>();

      if (node["experiment-label-paper"])
        m_config.experiment_label_paper = node["experiment-label-paper"].as<bool>();

      if (node["extra-label"])
        m_config.extra_label = node["extra-label"].as<std::string>();

      if (node["luminosity-label"])
        m_config.lumi_label = node["luminosity-label"].as<std::string>();

      if (node["root"]) {
        // Use # to remove root path - for FR2 combination from different folders
        if (node["root"].as<std::string>() == "#")
            m_config.root = "";
        else
            m_config.root = fs::absolute(fs::path(node["root"].as<std::string>()), histogramsPath).string();
      }

      if (node["scale"])
        m_config.scale = node["scale"].as<float>();

      if (node["eras"]) {
        auto availEras = node["eras"].as<std::vector<std::string>>();
        if ( CommandLineCfg::get().era.empty() ) {
          m_config.eras = availEras;
        } else {
          const auto reqEra = CommandLineCfg::get().era;
          if ( std::end(availEras) == std::find(std::begin(availEras), std::end(availEras), reqEra) ) {
            throw std::runtime_error("Requested era "+reqEra+" not found in configuration file");
          }
          m_config.eras = { CommandLineCfg::get().era };
        }
      }

      if (node["luminosity"]) {
        const auto& lumiNd = node["luminosity"];
        if ( lumiNd.IsScalar() ) {
          m_config.luminosity[""] = lumiNd.as<float>();
        } else if ( lumiNd.IsMap() ) {
          float totLumi = 0;
          for ( const auto& era : m_config.eras ) {
            const auto eraLumi = lumiNd[era].as<float>();
            m_config.luminosity[era] = eraLumi;
            totLumi += eraLumi;
          }
          m_config.luminosity[""] = totLumi;
        } else {
          throw YAML::ParserException(YAML::Mark::null_mark(), "luminosity should be a single value or a map (one value per era)");
        }
      } else {
        throw YAML::ParserException(YAML::Mark::null_mark(), "'configuration' block is missing luminosity");
      }

      if (node["no-lumi-rescaling"])
        m_config.no_lumi_rescaling = node["no-lumi-rescaling"].as<bool>();

      if (node["luminosity-error"]) {
        float value = node["luminosity-error"].as<float>();

        if (value > 0) {
          // Create a 'luminosity' systematic error
          YAML::Node syst;
          syst["type"] = "const";
          syst["pretty-name"] = "Luminosity";
          syst["value"] = value + 1;

          YAML::Node syst_node;
          syst_node["lumi"] = syst;

          f["systematics"].push_back(syst_node);
        }
      }

      if (node["error-fill-color"])
        m_config.error_fill_color = loadColor(node["error-fill-color"]);

      if (node["error-fill-style"])
        m_config.error_fill_style = loadColor(node["error-fill-style"]);

      if (node["staterror-fill-color"])
        m_config.staterror_fill_color = loadColor(node["staterror-fill-color"]);

      if (node["staterror-fill-style"])
        m_config.staterror_fill_style = loadColor(node["staterror-fill-style"]);

      if (node["fit-line-style"])
        m_config.fit_line_style = node["fit-line-style"].as<int16_t>();

      if (node["fit-line-width"])
        m_config.fit_line_width = node["fit-line-width"].as<int16_t>();

      if (node["fit-line-color"])
        m_config.fit_line_color = loadColor(node["fit-line-color"]);

      if (node["fit-error-fill-style"])
        m_config.fit_error_fill_style = node["fit-error-fill-style"].as<int16_t>();

      if (node["fit-error-fill-color"])
        m_config.fit_error_fill_color = loadColor(node["fit-error-fill-color"]);

      if (node["fit-n-points"])
        m_config.fit_n_points = node["fit-n-points"].as<uint16_t>();

      if (node["ratio-fit-line-style"])
        m_config.ratio_fit_line_style = node["ratio-fit-line-style"].as<int16_t>();

      if (node["ratio-fit-line-width"])
        m_config.ratio_fit_line_width = node["ratio-fit-line-width"].as<int16_t>();

      if (node["ratio-fit-line-color"])
        m_config.ratio_fit_line_color = loadColor(node["ratio-fit-line-color"]);

      if (node["ratio-fit-error-fill-style"])
        m_config.ratio_fit_error_fill_style = node["ratio-fit-error-fill-style"].as<int16_t>();

      if (node["ratio-fit-error-fill-color"])
        m_config.ratio_fit_error_fill_color = loadColor(node["ratio-fit-error-fill-color"]);

      if (node["ratio-fit-n-points"])
        m_config.ratio_fit_n_points = node["ratio-fit-n-points"].as<uint16_t>();

      if (node["blinded-range-fill-color"])
        m_config.blinded_range_fill_color = loadColor(node["blinded-range-fill-color"]);

      if (node["blinded-range-fill-style"])
        m_config.blinded_range_fill_style = node["blinded-range-fill-style"].as<uint16_t>();

      if (node["uncertainty-label"])
        m_config.uncertainty_label = node["uncertainty-label"].as<std::string>();

      if (node["syst-only"])
        m_config.syst_only = node["syst-only"].as<bool>();

      m_config.line_style.parse(node);

      if (node["labels"]) {
        YAML::Node labels = node["labels"];
        m_config.labels = parseLabelsNode(labels);
      }

      if (node["y-axis-format"])
        m_config.y_axis_format = node["y-axis-format"].as<std::string>();

      if (node["ratio-y-axis"])
        m_config.ratio_y_axis_title = node["ratio-y-axis"].as<std::string>();

      if (node["ratio-style"])
        m_config.ratio_style = node["ratio-style"].as<std::string>();

      if (node["mode"])
          m_config.mode = node["mode"].as<std::string>();

      if (node["tree-name"])
          m_config.tree_name = node["tree-name"].as<std::string>();

      if (node["transparent-background"])
          m_config.transparent_background = node["transparent-background"].as<bool>();

      if (node["show-overflow"])
          m_config.show_overflow = node["show-overflow"].as<bool>();

      if (node["show-onlyoverflow"]) {
          m_config.show_onlyoverflow = node["show-onlyoverflow"].as<bool>();
          m_config.show_overflow = False;
      }

      if (node["errors-type"])
          m_config.errors_type = string_to_errors_type(node["errors-type"].as<std::string>());

      if (node["yields-table-stretch"])
        m_config.yields_table_stretch = node["yields-table-stretch"].as<float>();

      if (node["yields-table-align"])
        m_config.yields_table_align = node["yields-table-align"].as<std::string>();

      if (node["yields-table-text-align"])
        m_config.yields_table_text_align = node["yields-table-text-align"].as<std::string>();

      if (node["yields-table-numerical-precision-yields"])
        m_config.yields_table_num_prec_yields = node["yields-table-numerical-precision-yields"].as<int>();

      if (node["yields-table-numerical-precision-ratio"])
        m_config.yields_table_num_prec_ratio = node["yields-table-numerical-precision-ratio"].as<int>();

      if (node["book-keeping-file"])
        m_config.book_keeping_file_name = node["book-keeping-file"].as<std::string>();

      // Axis size
      if (node["x-axis-label-size"])
        m_config.x_axis_label_size = node["x-axis-label-size"].as<float>();

      if (node["y-axis-label-size"])
        m_config.y_axis_label_size = node["y-axis-label-size"].as<float>();

      if (node["x-axis-top-ticks"])
        m_config.x_axis_top_ticks = node["x-axis-top-ticks"].as<bool>();

      if (node["y-axis-right-ticks"])
        m_config.y_axis_right_ticks = node["y-axis-right-ticks"].as<bool>();

      if (node["generated-events-histogram"])
        m_config.generated_events_histogram = node["generated-events-histogram"].as<std::string>();

      if (node["generated-events-bin"])
        m_config.generated_events_bin = node["generated-events-bin"].as<int>();
    }

    // Retrieve files/processes configuration
    YAML::Node files = f["files"];

    size_t process_id = 0;
    for (YAML::const_iterator it = files.begin(); it != files.end(); ++it) {
        File file;
        if (!CommandLineCfg::get().do_qcd and (it->first.as<std::string>()).find("QCD") != std::string::npos) continue;
        //Custom switches
        if (!CommandLineCfg::get().dyincl and (it->first.as<std::string>()).find("DYJetsToLL_M50_amc") != std::string::npos) continue;
        else if (CommandLineCfg::get().dyincl and (it->first.as<std::string>()).find("DYJetsToLL_M50_HT") != std::string::npos) continue;

        std::string fname_ = it->first.as<std::string>();
        std::transform(fname_.begin(), fname_.end(), fname_.begin(), [](unsigned char c){ return std::tolower(c); });
        bool isLfvSignal = fname_.find("_lfv_") != std::string::npos;

        if(!CommandLineCfg::get().allSig and isLfvSignal and (it->first.as<std::string>()).find("Vector") == std::string::npos ) continue;
        // This option must be used with allSig
        if(CommandLineCfg::get().selectSig.length() > 0 and isLfvSignal and (it->first.as<std::string>()).find(CommandLineCfg::get().selectSig) == std::string::npos ) continue;

        if(CommandLineCfg::get().noSig and isLfvSignal) continue;

        if(CommandLineCfg::get().binned and (((it->first.as<std::string>()).find("ttbarsignal") != std::string::npos
                or (it->first.as<std::string>()).find("ttbartwbb4l") != std::string::npos
                or (it->first.as<std::string>()).find("singletw") != std::string::npos)
            and (it->first.as<std::string>()).find("Bin") == std::string::npos)) continue;
        else if(!CommandLineCfg::get().binned
            and (((it->first.as<std::string>()).find("ttbarsignal") != std::string::npos
                or (it->first.as<std::string>()).find("ttbartwbb4l") != std::string::npos
                or (it->first.as<std::string>()).find("singletw") != std::string::npos)
            and (it->first.as<std::string>()).find("Bin") != std::string::npos)) continue;

        if (files.Type() == YAML::NodeType::Map)
            parseFileNode(file, it->first, it->second);
        else
            parseFileNode(file, *it);

        file.id = process_id++;
        if ( filter_eras(file) ) {
          m_files.push_back(file);
        }
    }

    if (! expandFiles())
        return false;

    std::sort(m_files.begin(), m_files.end(), [](const File& a, const File& b) {
      return a.order < b.order;
     });

    YAML::Node legend_groups = f["groups"];

    for (YAML::const_iterator it = legend_groups.begin(); it != legend_groups.end(); ++it) {
      Group group;

      group.name = it->first.as<std::string>();

      YAML::Node node = it->second;

      // Find the first file belonging to this group, and use its type to set
      // default style values
      const auto file = std::find_if(m_files.begin(), m_files.end(), [&group](const File& file) {
          return file.legend_group == group.name;
        });

      // Is this group actually used?
      if (file == m_files.end())
          continue;

      group.plot_style = std::make_shared<PlotStyle>();
      group.plot_style->loadFromYAML(node, file->type);

      m_legend_groups[group.name] = group;

      if ( node["order"] ) {
        const auto groupOrder = node["order"].as<int16_t>();
        for ( auto& file : m_files ) {
          if ( ( file.legend_group == group.name ) && ( file.order == std::numeric_limits<int16_t>::min() ) ) {
            file.order = groupOrder;
          }
        }
      }
    }

    std::sort(m_files.begin(), m_files.end(), [](const File& a, const File& b) {
      return a.order < b.order;
     });

    // Remove non-existant groups from files and update yields group
    for (auto& file: m_files) {
      if (!file.legend_group.empty() && !m_legend_groups.count(file.legend_group)) {
        std::cout << "Warning: group " << file.legend_group << " (used for file " << file.pretty_name << ") not found, ignoring" << std::endl;
        file.legend_group = "";
      }

      if (file.yields_group.empty()) {
        if (!file.legend_group.empty()) {
            file.yields_group = m_legend_groups[file.legend_group].plot_style->legend;
        }

        if (file.yields_group.empty())
            file.yields_group = file.plot_style->legend;

        if (file.yields_group.empty())
            file.yields_group = file.pretty_name;
      }
    }

    // List systematics
    if (f["systematics"]) {
        YAML::Node systs = f["systematics"];

        for (YAML::const_iterator it = systs.begin(); it != systs.end(); ++it) {
            parseSystematicsNode(*it);
        }
    }

    // Retrieve plots configuration
    if (! f["plots"]) {
      throw YAML::ParserException(YAML::Mark::null_mark(), "You must specify at least one plot in your configuration file");
    }

    YAML::Node plots = f["plots"];

    for (YAML::const_iterator it = plots.begin(); it != plots.end(); ++it) {
      Plot plot;

      plot.name = it->first.as<std::string>();

      YAML::Node node = it->second;
      if (node["exclude"])
        plot.exclude = node["exclude"].as<std::string>();

      if (node["x-axis"])
        plot.x_axis = node["x-axis"].as<std::string>();

      if (node["y-axis"])
        plot.y_axis = node["y-axis"].as<std::string>();

      if (node["ratio-y-axis"])
        plot.ratio_y_axis_title = node["ratio-y-axis"].as<std::string>();
      else
        plot.ratio_y_axis_title = m_config.ratio_y_axis_title;

      plot.y_axis_format = m_config.y_axis_format;
      if (node["y-axis-format"])
        plot.y_axis_format = node["y-axis-format"].as<std::string>();

      if (node["normalized"])
        plot.normalized = node["normalized"].as<bool>();

      if (node["signal-normalize-data"])
        plot.signal_normalize_data = node["signal-normalize-data"].as<bool>();

      if (node["no-data"])
        plot.no_data = node["no-data"].as<bool>();

      if (node["override"])
        plot.override = node["override"].as<bool>();

      Log log_y = False;
      if (node["log-y"]) {
        log_y = parse_log(node["log-y"]);
      }
      if (log_y != Both)
        plot.log_y = (bool) log_y;

      Log log_x = False;
      if (node["log-x"]) {
        log_x = parse_log(node["log-x"]);
      }
      if (log_x != Both)
        plot.log_x = (bool) log_x;

      if (node["save-extensions"])
        plot.save_extensions = node["save-extensions"].as<std::vector<std::string>>();

      if (node["show-ratio"])
        plot.show_ratio = node["show-ratio"].as<bool>();

      if (node["fit-ratio"])
        plot.fit_ratio = node["fit-ratio"].as<bool>();

      if (node["fit"])
        plot.fit = node["fit"].as<bool>();

      if (node["fit-function"])
        plot.fit_function = node["fit-function"].as<std::string>();

      if (node["fit-legend"])
        plot.fit_legend = node["fit-legend"].as<std::string>();

      if (node["fit-legend-position"])
        plot.fit_legend_position = node["fit-legend-position"].as<Point>();

      if (node["fit-range"])
        plot.fit_range = node["fit-range"].as<Range>();

      if (node["ratio-fit-function"])
        plot.ratio_fit_function = node["ratio-fit-function"].as<std::string>();

      if (node["ratio-fit-legend"])
        plot.ratio_fit_legend = node["ratio-fit-legend"].as<std::string>();

      if (node["ratio-fit-legend-position"])
        plot.ratio_fit_legend_position = node["ratio-fit-legend-position"].as<Point>();

      if (node["ratio-fit-range"])
        plot.ratio_fit_range = node["ratio-fit-range"].as<Range>();

      if (node["show-errors"])
        plot.show_errors = node["show-errors"].as<bool>();

      if (node["x-axis-range"])
        plot.x_axis_range = node["x-axis-range"].as<Range>();
      plot.log_x_axis_range = plot.x_axis_range;

      if (node["log-x-axis-range"])
        plot.log_x_axis_range = node["log-x-axis-range"].as<Range>();

      if (node["y-axis-auto-range"])
        plot.y_axis_auto_range = node["y-axis-auto-range"].as<bool>();

      if (node["y-axis-range"])
        plot.y_axis_range = node["y-axis-range"].as<Range>();
      plot.log_y_axis_range = plot.y_axis_range;

      if (node["log-y-axis-range"])
        plot.log_y_axis_range = node["log-y-axis-range"].as<Range>();

      if (node["ratio-y-axis-range"])
        plot.ratio_y_axis_range = node["ratio-y-axis-range"].as<Range>();

      if (node["ratio-y-axis-auto-range"])
        plot.ratio_y_axis_auto_range = node["ratio-y-axis-auto-range"].as<bool>();

      if (node["ratio-draw-mcstat-error"])
        plot.ratio_draw_mcstat_error = node["ratio-draw-mcstat-error"].as<bool>();

      if (node["draw-siglike-unc"])
        plot.draw_siglike_unc = node["draw-siglike-unc"].as<bool>();

      if (node["post-fit"])
        plot.post_fit = node["post-fit"].as<bool>();

      if (node["blinded-range"])
        plot.blinded_range = node["blinded-range"].as<Range>();

      if (node["y-axis-show-zero"])
        plot.y_axis_show_zero = node["y-axis-show-zero"].as<bool>();

      if (node["inherits-from"])
        plot.inherits_from = node["inherits-from"].as<std::string>();

      if (node["rebin"])
        plot.rebin = node["rebin"].as<uint16_t>();

      if (node["rebin-arr"])
        plot.rebin_arr = node["rebin-arr"].as<std::vector<double>>();

      if (node["labels"]) {
        YAML::Node labels = node["labels"];
        plot.labels = parseLabelsNode(labels);
      }

      // Change legend by hand, currently can change only one entry FIXME
      if (node["change-legend"]) {
        plot.change_legend = node["change-legend"].as<bool>();
        plot.legend_name_org = node["legend-name-org"].as<std::string>();
        plot.legend_name_new = node["legend-name-new"].as<std::string>();
      }

      if (node["extra-label"])
        plot.extra_label = node["extra-label"].as<std::string>();

      if (node["legend-position"])
        plot.legend_position = node["legend-position"].as<Position>();
      else
        plot.legend_position = m_legend.position;

      if (node["legend-columns"])
        plot.legend_columns = node["legend-columns"].as<size_t>();
      else
        plot.legend_columns = m_legend.columns;

      if (node["show-overflow"])
        plot.show_overflow = node["show-overflow"].as<bool>();
      else if (node["show-onlyoverflow"]) {
        plot.show_onlyoverflow = node["show-onlyoverflow"].as<bool>();
        plot.show_overflow = False;
      } else
        plot.show_overflow = m_config.show_overflow;

      if (node["errors-type"])
        plot.errors_type = string_to_errors_type(node["errors-type"].as<std::string>());
      else
        plot.errors_type = m_config.errors_type;

      if (node["binning-x"])
        plot.binning_x = node["binning-x"].as<uint16_t>();

      if (node["binning-y"])
        plot.binning_y = node["binning-y"].as<uint16_t>();

      if (node["draw-string"])
        plot.draw_string = node["draw-string"].as<std::string>();

      if (node["selection-string"])
        plot.selection_string = node["selection-string"].as<std::string>();

      if (node["for-yields"])
        plot.use_for_yields = node["for-yields"].as<bool>();

      if (node["yields-title"])
        plot.yields_title = node["yields-title"].as<std::string>();
      else
        plot.yields_title = plot.name;

      if (node["yields-table-order"])
        plot.yields_table_order = node["yields-table-order"].as<int>();

      if (node["vertical-lines"]) {
        for (const auto& line: node["vertical-lines"]) {
          plot.lines.push_back(Line(line, VERTICAL));
        }
      }

      if (node["horizontal-lines"]) {
        for (const auto& line: node["horizontal-lines"]) {
          plot.lines.push_back(Line(line, HORIZONTAL));
        }
      }

      if (node["lines"]) {
        for (const auto& line: node["lines"]) {
          plot.lines.push_back(Line(line, UNSPECIFIED));
        }
      }

      for (auto& line: plot.lines) {
        if (! line.style)
          line.style = m_config.line_style;
      }

      if (node["book-keeping-folder"]) {
        plot.book_keeping_folder = node["book-keeping-folder"].as<std::string>();
      }

      plot.renaming_ops = parseRenameNode(node);

      if (node["sort-by-yields"]) {
        plot.sort_by_yields = node["sort-by-yields"].as<bool>();
      }

      // Axis size
      if (node["x-axis-label-size"])
        plot.x_axis_label_size = node["x-axis-label-size"].as<float>();
      else
        plot.x_axis_label_size = m_config.x_axis_label_size;

      if (node["y-axis-label-size"])
        plot.y_axis_label_size = node["y-axis-label-size"].as<float>();
      else
        plot.y_axis_label_size = m_config.y_axis_label_size;

      // Show or hide ticks
      if (node["x-axis-hide-ticks"])
        plot.x_axis_hide_ticks = node["x-axis-hide-ticks"].as<bool>();

      if (node["y-axis-hide-ticks"])
        plot.y_axis_hide_ticks = node["y-axis-hide-ticks"].as<bool>();

      // Additional option for plot scaling
      if (node["scale-option"])
        plot.scale_option = node["scale-option"].as<std::string>();

      // Handle log
      std::vector<bool> logs_x;
      std::vector<bool> logs_y;

      if (log_x == Both) {
        logs_x.insert(logs_x.end(), {false, true});
      } else {
        logs_x.push_back(plot.log_x);
      }

      if (log_y == Both) {
        logs_y.insert(logs_y.end(), {false, true});
      } else {
        logs_y.push_back(plot.log_y);
      }

      int log_counter(0);
      for (auto x: logs_x) {
        for (auto y: logs_y) {
          Plot p = plot;
          p.log_x = x;
          p.log_y = y;
          // If the plot is used for yields, they should be output only once
          if(log_counter && plot.use_for_yields)
            p.use_for_yields = false;

          if (p.log_x)
            p.output_suffix += "_logx";

          if (p.log_y)
            p.output_suffix += "_logy";

          m_plots.push_back(p);
          ++log_counter;
        }
      }
    }

    // If at least one plot has 'override' set to true, keep only plots which do
    if( std::find_if(m_plots.begin(), m_plots.end(), [](Plot &plot){ return plot.override; }) != m_plots.end() ){
      auto new_end = std::remove_if(m_plots.begin(), m_plots.end(), [](Plot &plot){ return !plot.override; });
      m_plots.erase(new_end, m_plots.end());
    }

    parseLumiLabel();

    if (CommandLineCfg::get().verbose) {
        std::cout << " done." << std::endl;
    }

    return true;
  }

  void plotIt::parseLumiLabel() {

    boost::format formatter = get_formatter(m_config.lumi_label);

    float lumi = m_config.luminosity[""] / 1000.;
    formatter % lumi;

    m_config.lumi_label = formatter.str();
  }

  void plotIt::fillLegend(TLegend& legend, const Plot& plot, bool with_uncertainties) {
      std::vector<LegendEntry> legend_entries[plot.legend_columns];

      auto getLegendEntryFromFile = [&](File& file, LegendEntry& entry) {
          if (file.legend_group.length() > 0 && m_legend_groups.count(file.legend_group) && m_legend_groups[file.legend_group].plot_style->legend.length() > 0) {
              if (m_legend_groups[file.legend_group].added)
                  return false;
              m_legend_groups[file.legend_group].added = true;

              const auto& plot_style = m_legend_groups[file.legend_group].plot_style;
              int16_t order = 0;
              if (plot_style->legend_order == 0) order = file.order;
              //entry = {file.object, plot_style->legend, plot_style->legend_style, plot_style->legend_order};
              entry = {file.object, plot_style->legend, plot_style->legend_style, order};
          } else if (file.plot_style.get() && file.plot_style->legend.length() > 0) {
              entry = {file.object, file.plot_style->legend, file.plot_style->legend_style, file.plot_style->legend_order};
          }

          return true;
      };

      auto getEntries = [&](Type type) {
          std::vector<LegendEntry> entries;
          for (File& file: m_files) {
              if (file.type == type) {
                  LegendEntry entry;
                  if (getLegendEntryFromFile(file, entry)) {
                      entries.push_back(entry);
                  }
              }
          }

          for (const auto& entry: m_config.static_legend_entries[type]) {
            entries.push_back(entry);
          }

          if (CommandLineCfg::get().desytop)
              std::sort(entries.begin(), entries.end(), [](const LegendEntry& a, const LegendEntry& b) { return a.order > b.order; });
          else
              std::sort(entries.begin(), entries.end(), [](const LegendEntry& a, const LegendEntry& b) { return a.order < b.order; });

          return entries;
      };
/*
      // First, add data, always on first column
      if (!plot.no_data) {
          std::vector<LegendEntry> entries = getEntries(DATA);
          for (const auto& entry: entries)
              legend_entries[0].push_back(entry);
      }

      // Then MC, spanning on the remaining columns
      size_t index = 0;
      std::vector<LegendEntry> entries = getEntries(MC);
      for (const LegendEntry& entry: entries) {
          size_t column_index = (plot.legend_columns == 1) ? 0 : ((index % (plot.legend_columns - 1)) + 1);
          legend_entries[column_index].push_back(entry);
          index++;
      }

      // Signal, also on the first column
      entries = getEntries(SIGNAL);
      for (const LegendEntry& entry: entries) {
          legend_entries[0].push_back(entry);
      }

      // Finally, if requested, the uncertainties entry
      if (with_uncertainties)
          legend_entries[0].push_back({m_config.uncertainty_label, "f", m_config.error_fill_style, m_config.error_fill_color, 0});
*/
      // First MC, spanning on the columns
      size_t index = 0;
      std::vector<LegendEntry> entries = getEntries(MC);
      for (const LegendEntry& entry: entries) {
          size_t column_index = index % (plot.legend_columns);
          legend_entries[column_index].push_back(entry);
          index++;
      }

      // Next, signals
      entries = getEntries(SIGNAL);
      for (const LegendEntry& entry: entries) {
          size_t column_index = index % (plot.legend_columns);
          legend_entries[column_index].push_back(entry);
          index++;
      }

      // Next, data
      if (!plot.no_data) {
          std::vector<LegendEntry> entries = getEntries(DATA);
          for (const auto& entry: entries){
              size_t column_index = index % (plot.legend_columns);
              legend_entries[column_index].push_back(entry);
              index++;
          }
      }

      // Finally, if requested, the uncertainties entry
      if (with_uncertainties){
          size_t column_index = index % (plot.legend_columns);
          legend_entries[column_index].push_back({m_config.uncertainty_label, "f", m_config.error_fill_style, m_config.error_fill_color, 0});
      }

      // Ensure all columns have the same size
      size_t max_size = 0;
      for (size_t i = 0; i < plot.legend_columns; i++) {
          max_size = std::max(max_size, legend_entries[i].size());
      }

      for (size_t i = 0; i < plot.legend_columns; i++) {
          legend_entries[i].resize(max_size, LegendEntry());
      }

      // Add entries to the legend
      for (size_t i = 0; i < (plot.legend_columns * max_size); i++) {
          size_t column_index = (i % plot.legend_columns);
          size_t row_index = static_cast<size_t>(i / static_cast<float>(plot.legend_columns));
          LegendEntry& entry = legend_entries[column_index][row_index];
          TLegendEntry* e = legend.AddEntry(entry.object, entry.legend.c_str(), entry.style.c_str());
          entry.stylize(e);
      }
  }

  bool plotIt::plot(Plot& plot) {
    std::cout << "Plotting '" << plot.name << "'" << std::endl;

    bool hasMC = false;
    bool hasData = false;
    bool hasSignal = false;
    bool hasLegend = false;
    // Open all files, and find histogram in each
    for (File& file: m_files) {
      if (! loadObject(file, plot)) {
        return false;
      }

      hasLegend |= getPlotStyle(file)->legend.length() > 0;
      hasData |= file.type == DATA;
      hasMC |= file.type == MC;
      hasSignal |= file.type == SIGNAL;
    }

    // Can contains '/' if the plot is inside a folder
    fs::path plot_path = plot.name + plot.output_suffix;
    std::string plot_name = plot_path.filename().string();

    // Create canvas
    TCanvas c(plot_name.c_str(), plot_name.c_str(), m_config.width, m_config.height);

    if (m_config.transparent_background) {
        c.SetFillStyle(4000);
        c.SetFrameFillStyle(4000);
    }

    if ( m_files.empty() ) {
      std::cout << "No files selected" << std::endl;
      return false;
    }

    boost::optional<Summary> summary = ::plotIt::plot(m_files[0], c, plot);

    if (! summary)
      return false;

    if (CommandLineCfg::get().verbose) {
      ConsoleSummaryPrinter printer;
      printer.print(*summary);
    }

    if (plot.log_y)
      c.SetLogy();

    //if (plot.log_x)
    //  c.SetLogx();

    Position legend_position = plot.legend_position;

    // Build legend
    TLegend legend(legend_position.x1, legend_position.y1, legend_position.x2, legend_position.y2);
    //legend.SetTextFont(62);
    legend.SetTextFont(42);
    //legend.SetTextSize(0.045);
    legend.SetTextSize(0.055);
    if (CommandLineCfg::get().desytop)
        legend.SetTextSize(0.045);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetNColumns(plot.legend_columns);

    fillLegend(legend, plot, hasMC && plot.show_errors);

    if (plot.change_legend) {
      TList *p = legend.GetListOfPrimitives();
      TIter next(p);
      TObject *obj;
      TLegendEntry *le;
      while ((obj = next())) {
        le = (TLegendEntry*)obj;
        std::string le_name = le->GetLabel();
        if (le_name == plot.legend_name_org) le->SetLabel((plot.legend_name_new).c_str());
      }
    }

    legend.Draw();

    float topMargin = m_config.margin_top;
    if (plot.show_ratio)
      topMargin /= .6666;

    // Move exponent label if shown
    TGaxis::SetMaxDigits(3);
    TGaxis::SetExponentOffset(-0.03, 0.01, "y");

    // Luminosity label
    if (m_config.lumi_label.length() > 0) {
      std::shared_ptr<TPaveText> pt = std::make_shared<TPaveText>(m_config.margin_left, 1 - 0.2 * topMargin, 1 - m_config.margin_right, 1, "NDC");
      TemporaryPool::get().add(pt);

      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->SetMargin(0);
      pt->SetTextFont(62);
      pt->SetTextSize(0.8 * topMargin);
      pt->SetTextAlign(33);

      pt->AddText(m_config.lumi_label.c_str());
      pt->Draw();
    }

    // Experiment
    if (m_config.experiment.length() > 0) {
      std::shared_ptr<TPaveText> pt;
      if (m_config.experiment_label_paper)
          //pt = std::make_shared<TPaveText>(1.15 * m_config.margin_left, 1 - 2.75 * topMargin, 1 - m_config.margin_right, 1, "brNDC");
          pt = std::make_shared<TPaveText>(1.15 * m_config.margin_left, 1 - 2.6 * topMargin, 1 - m_config.margin_right, 1 - 1.0 * topMargin, "brNDC");
      else pt = std::make_shared<TPaveText>(m_config.margin_left, 1 - 0.5 * topMargin, 1 - m_config.margin_right, 1, "brNDC");
      TemporaryPool::get().add(pt);

      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->SetMargin(0);
      pt->SetTextFont(62);
      if (m_config.experiment_label_paper)
          pt->SetTextSize(0.9 * topMargin);
      else pt->SetTextSize(0.65 * topMargin);
      pt->SetTextAlign(13);

      std::string text = m_config.experiment;
      std::string text2 = m_config.extra_label;
      if (m_config.extra_label.length() || plot.extra_label.length()) {
        std::string extra_label = plot.extra_label;
        if (extra_label.length() == 0) {
          extra_label = m_config.extra_label;
        }

        //boost::format fmt("%s #font[52]{#scale[0.76]{%s}}");
        //fmt % m_config.experiment % extra_label;

        //text = fmt.str();

        boost::format fmt("#font[52]{#scale[0.58]{%s}}");
        fmt % extra_label;

        text2 = fmt.str();
      }

      pt->AddText(text.c_str());
      pt->AddText(text2.c_str());
      pt->Draw();
    }

    c.cd();

    const auto& labels = mergeLabels(plot.labels);

    // Labels
    for (auto& label: labels) {

      std::shared_ptr<TLatex> t(new TLatex(label.position.x, label.position.y, label.text.c_str()));
      t->SetNDC(true);
      //t->SetTextFont(64);
      t->SetTextFont(label.font);
      t->SetTextSize(label.size);
      t->Draw();

      TemporaryPool::get().add(t);
    }

    fs::path rootDir = m_outputPath;
    fs::path outputName = rootDir / plot_path;

    // Ensure path exists
    fs::create_directories(outputName.parent_path());

    for (const std::string& extension: plot.save_extensions) {
      fs::path plotPathWithExtension = plot_path.replace_extension(extension);

      std::string finalPlotPathWithExtension = applyRenaming(plot.renaming_ops, plotPathWithExtension.native());
      fs::path finalOutputName = rootDir / finalPlotPathWithExtension;

      c.SaveAs(finalOutputName.c_str());
    }

    if (m_config.book_keeping_file) {
      TDirectory* root = m_config.book_keeping_file.get();
      if (!plot.book_keeping_folder.empty() || !plot_path.parent_path().empty()) {
        // Look in the cache if we have this folder. This avoid querying the file each time we save a plot
        std::string path = (!plot.book_keeping_folder.empty()) ? plot.book_keeping_folder : plot_path.parent_path().string();
        auto it = m_book_keeping_folders.find(path);
        if (it == m_book_keeping_folders.end()) {
          root = ::plotIt::getDirectory(m_config.book_keeping_file.get(), path);
          m_book_keeping_folders.emplace(path, root);
        } else {
          root = it->second;
        }
      }
      root->WriteTObject(&c, nullptr, "Overwrite");
    }

    // Clean all temporary resources
    TemporaryPool::get().clear();

    // Reset groups
    for (auto& group: m_legend_groups) {
      group.second.added = false;
    }

    return true;
  }


  // yield table
  bool plotIt::yields(std::vector<Plot>::iterator plots_begin, std::vector<Plot>::iterator plots_end){
    std::cout << "Producing LaTeX yield table.\n";

    std::map<std::string, double> data_yields;

    std::map< std::string, std::map<std::string, std::pair<double, double> > > mc_yields;
    std::map< std::string, double > mc_total;
    std::map< std::string, double > mc_total_sqerrs;
    std::set<std::string> mc_processes;

    std::map< std::string, std::map<std::string, std::pair<double, double> > > signal_yields;
    std::set<std::string> signal_processes;

    std::map<
        std::tuple<Type, std::string, std::string>, // Type, category, systematics name
        double
    > process_systematics;

    std::map<
        std::string,
        std::map<Type, double>
    > total_systematics_squared;

    std::vector< std::pair<int, std::string> > categories;

    bool has_data(false);

    for ( auto it = plots_begin; it != plots_end; ++it ) {
      auto& plot = *it;
      if (!plot.use_for_yields)
        continue;

      if (plot.yields_title.find("$") == std::string::npos)
          replace_substr(plot.yields_title, "_", "\\_");

      if( std::find_if(categories.begin(), categories.end(), [&](const std::pair<int, std::string> &x){ return x.second == plot.yields_title; }) != categories.end() )
          continue;
      categories.push_back( std::make_pair(plot.yields_table_order, plot.yields_title) );

      std::map<std::tuple<Type, std::string>, double> plot_total_systematics;

      // Open all files, and find histogram in each
      for (File& file: m_files) {
        if (! loadObject(file, plot)) {
          std::cout << "Could not retrieve plot from " << file.path << std::endl;
          return false;
        }

        if ( file.type == DATA ){
          TH1* h = dynamic_cast<TH1*>(file.object);
          data_yields[plot.yields_title] += h->Integral(0, h->GetNbinsX() + 1);
          has_data = true;
          continue;
        }

        std::string process_name = file.yields_group;

        if (process_name.find("$") == std::string::npos)
            replace_substr(process_name, "_", "\\_");

        if (process_name.find("#") != std::string::npos) {
            // We assume it's a ROOT LaTeX string. Enclose the string into $$, and replace
            // '#' by '\'

            replace_substr(process_name, "#", R"(\)");
            process_name = "$" + process_name + "$";
        }

        std::pair<double, double> yield_sqerror;
        TH1* hist( dynamic_cast<TH1*>(file.object) );

        if (m_config.generated_events_histogram.length() > 0 and file.generated_events < 0 and file.type != DATA ) {
            //generated_events = - 1.0 if not declared in file yaml
            std::shared_ptr<TFile> input(TFile::Open(file.path.c_str()));
            TH1* hevt = dynamic_cast<TH1*>(input->Get(m_config.generated_events_histogram.c_str()));
            file.generated_events = hevt->GetBinContent(m_config.generated_events_bin);
        }
        double factor = file.cross_section * file.branching_ratio / file.generated_events;

        if (! m_config.no_lumi_rescaling) {
          factor *= m_config.luminosity.at(file.era);
        }
        if (!CommandLineCfg::get().ignore_scales)
          factor *= m_config.scale * file.scale;

        if (!plot.is_rescaled)
          hist->Scale(factor);

        for (auto& syst: *file.systematics) {
          syst.update();
          syst.scale(factor);
        }

        // Retrieve yield and stat. error, taking overflow into account
        yield_sqerror.first = hist->IntegralAndError(0, hist->GetNbinsX() + 1, yield_sqerror.second);
        yield_sqerror.second = std::pow(yield_sqerror.second, 2);

        // Add systematics
        double file_total_systematics = 0;
        for (auto& syst: *file.systematics) {

          TH1* nominal_shape = static_cast<TH1*>(syst.nominal_shape.get());
          TH1* up_shape = static_cast<TH1*>(syst.up_shape.get());
          TH1* down_shape = static_cast<TH1*>(syst.down_shape.get());

          if (! nominal_shape || ! up_shape || ! down_shape)
              continue;

          double nominal_integral = nominal_shape->Integral(0, nominal_shape->GetNbinsX() + 1);
          double up_integral = up_shape->Integral(0, up_shape->GetNbinsX() + 1);
          double down_integral = down_shape->Integral(0, down_shape->GetNbinsX() + 1);

          double total_syst_error = std::max(
                  std::abs(up_integral - nominal_integral),
                  std::abs(nominal_integral - down_integral)
          );

          file_total_systematics += total_syst_error * total_syst_error;

          auto key = std::make_tuple(file.type, syst.name());
          plot_total_systematics[key] += total_syst_error;
        }

        // file_total_systematics contains the quadratic sum of all the systematics for this file
        process_systematics[std::make_tuple(file.type, plot.yields_title, process_name)] += std::sqrt(file_total_systematics);

        if ( file.type == MC ){
          ADD_PAIRS(mc_yields[plot.yields_title][process_name], yield_sqerror);
          mc_total[plot.yields_title] += yield_sqerror.first;
          mc_total_sqerrs[plot.yields_title] += yield_sqerror.second;
          mc_processes.emplace(process_name);
        }
        if ( file.type == SIGNAL ){
          ADD_PAIRS(signal_yields[plot.yields_title][process_name], yield_sqerror);
          signal_processes.emplace(process_name);
        }
      }

      // Get the total systematics for this category
      for (auto& syst: plot_total_systematics) {
        total_systematics_squared[plot.yields_title][std::get<0>(syst.first)] += syst.second * syst.second;
      }
    }

    if( ( !(mc_processes.size()+signal_processes.size()) && !has_data ) || !categories.size() ){
      std::cout << "No processes/data/categories defined\n";
      return false;
    }

    //std::vector<std::string> mc_processes_sorted;
    std::map< std::string, std::string > mc_processes_sorted;
    for(auto &process_name: mc_processes) {
      auto process_name_modified = process_name;
      if (process_name.find_first_of("0123456789") != std::string::npos) {
          // trick to sort process: put integer in front of the process name
          // then we remove here
          // Assume we have < 10 processes

          auto begin = process_name.find_first_of("0123456789");
          process_name_modified = process_name.substr(begin+1);

          std::string mathmark = "$";
          if (process_name.compare(0, 1, mathmark) == 0) process_name_modified = "$" + process_name_modified;
      }
      //mc_processes_sorted.emplace_back(process_name_modified);
      mc_processes_sorted[process_name] = process_name_modified;
    }

    std::map< std::string, std::string > signal_processes_sorted;
    for(auto &process_name: signal_processes) {
      auto process_name_modified = process_name;
      if (process_name.find_first_of("0123456789") != std::string::npos) {

          auto begin = process_name.find_first_of("0123456789");
          process_name_modified = process_name.substr(begin+1);

          std::string mathmark = "$";
          if (process_name.compare(0, 1, mathmark) == 0) process_name_modified = "$" + process_name_modified;
      }
      signal_processes_sorted[process_name] = process_name_modified;
    }

    // Sort according to user-defined order
    std::sort(categories.begin(), categories.end(), [](const std::pair<int, std::string>& cat1, const std::pair<int, std::string>& cat2){  return cat1.first < cat2.first; });

    std::ostringstream latexString;
    std::string tab("    ");

    latexString << std::setiosflags(std::ios_base::fixed);

    auto format_number_with_errors = [](double number, double error_low, double error_high, uint8_t number_precision, uint8_t error_precision) -> std::string {
        std::stringstream ss;
        ss << std::setiosflags(std::ios_base::fixed);
        ss << std::setprecision(number_precision);
        if (std::abs(error_high - error_low) > std::pow(10, -1 * error_precision)) {
            // Errors are really asymmetric
            ss << "$" << number << std::setprecision(error_precision) << "^{+" << error_high << "}_{-" << error_low << "}$";
        } else {
            // Symmetric errors
            //ss << "$" << number << R"( {\scriptstyle\ \pm\ )" << std::setprecision(error_precision) << error_low << "}$";
            ss << "$" << number << "$ &$\\pm$& $" << std::setprecision(error_precision) << error_low << "$";
        }

        return ss.str();
    };

    if( m_config.yields_table_align.find("h") != std::string::npos ){

      latexString << "\\renewcommand{\\arraystretch}{" << m_config.yields_table_stretch << "}\n";
      latexString << "\\begin{tabular}{ |l||";

      // tabular config.
      for(size_t i = 0; i < signal_processes.size(); ++i)
        latexString << m_config.yields_table_text_align << "|";
      if(signal_processes.size())
        latexString << "|";
      for(size_t i = 0; i < mc_processes.size(); ++i)
        latexString << m_config.yields_table_text_align << "|";
      if(mc_processes.size())
        latexString << "|" + m_config.yields_table_text_align << "||";
      if(has_data)
        latexString << m_config.yields_table_text_align << "||";
      if(has_data && mc_processes.size())
        latexString << m_config.yields_table_text_align << "||";
      latexString.seekp(latexString.tellp() - 2l);
      latexString << "| }\n" << tab << tab << "\\hline\n";

      // title line
      latexString << "    Cat. & ";
      for(auto &proc: signal_processes_sorted)
        latexString << proc.second << " & ";
      //for(auto &proc: mc_processes)
      for(auto &proc: mc_processes_sorted)
        latexString << proc.second << " & ";
      if( mc_processes.size() )
        latexString << "Tot. MC & ";
      if( has_data )
        latexString << "Data & ";
      if( has_data && mc_processes.size() )
        latexString << "Data/MC & ";
      latexString.seekp(latexString.tellp() - 2l);
      latexString << "\\\\\n" << tab << tab << "\\hline\n";

      // loop over each category
      for(auto& cat_pair: categories){

        std::string categ(cat_pair.second);
        latexString << tab << categ << " & ";
        latexString << std::setprecision(m_config.yields_table_num_prec_yields);

        for(auto &proc: signal_processes_sorted)
          latexString << "$" << signal_yields[categ][proc.first].first << " \\pm " << std::sqrt(signal_yields[categ][proc.first].second) << "$ & ";
          //latexString << "$" << signal_yields[categ][proc].first << " \\pm " << std::sqrt(signal_yields[categ][proc].second + std::pow(process_systematics[std::make_tuple(SIGNAL, categ, proc)], 2)) << "$ & ";

        //for(auto &proc: mc_processes) {
        for(auto &proc: mc_processes_sorted) {
          if(mc_yields[categ][proc.first].first > 0)
            latexString << "$" << mc_yields[categ][proc.first].first << " \\pm " << std::sqrt(mc_yields[categ][proc.first].second) << "$ & ";
            //latexString << "$" << mc_yields[categ][proc].first << " \\pm " << std::sqrt(mc_yields[categ][proc].second + std::pow(process_systematics[std::make_tuple(MC, categ, proc)], 2)) << "$ & ";
          else
            latexString << "$" << "0" << " \\pm " << std::sqrt(mc_yields[categ][proc.first].second) << "$ & ";
        }
        if( mc_processes.size() )
          latexString << "$" << mc_total[categ] << " \\pm " << std::sqrt(mc_total_sqerrs[categ]) << "$ & ";
          //latexString << "$" << mc_total[categ] << " \\pm " << std::sqrt(mc_total_sqerrs[categ] + total_systematics_squared[categ][MC]) << "$ & ";

        if( has_data ) {
          static const double alpha = 1. - 0.682689492;
          uint64_t yield = data_yields[cat_pair.second];
          double error_low = yield - ROOT::Math::gamma_quantile(alpha / 2., yield, 1.);
          double error_high = ROOT::Math::gamma_quantile_c(alpha / 2., yield, 1.) - yield;
          latexString << format_number_with_errors(yield, error_low, error_high, 0, m_config.yields_table_num_prec_yields) << " & ";
        }

        if( has_data && mc_processes.size() ){
          uint64_t data_yield = data_yields[categ];
          double ratio = data_yield / mc_total[categ];

          static const double alpha = 1. - 0.682689492;
          double error_data_low = data_yield - ROOT::Math::gamma_quantile(alpha / 2., data_yield, 1.);
          double error_data_high = ROOT::Math::gamma_quantile_c(alpha / 2., data_yield, 1.) - data_yield;

          double error_mc = std::sqrt(mc_total_sqerrs[categ] + total_systematics_squared[categ][MC]);

          double error_low = ratio * std::sqrt(std::pow(error_data_low / data_yields[categ], 2) +  std::pow(error_mc / mc_total[categ], 2));
          double error_high = ratio * std::sqrt(std::pow(error_data_high / data_yields[categ], 2) +  std::pow(error_mc / mc_total[categ], 2));

          //latexString << format_number_with_errors(ratio, error_low, error_high, m_config.yields_table_num_prec_ratio, m_config.yields_table_num_prec_ratio) << " & ";
          latexString << format_number_with_errors(ratio, error_low, error_high, m_config.yields_table_num_prec_ratio, 1) << " & ";
        }

        latexString.seekp(latexString.tellp() - 2l);
        latexString << "\\\\\n";
      }

      latexString << tab << tab << "\\hline\n\\end{tabular}\n";

    } else if (m_config.yields_table_align.find("v") != std::string::npos) {

        // Tabular header
        latexString << R"(\begin{tabular}{@{}lc)";
        std::string header = " & & ";
        for (size_t i = 0; i < categories.size(); i++) {
            latexString << "rcl";
            header += R"(\multicolumn{3}{c}{)" + categories[i].second + R"(})";
            if (i != (categories.size() - 1))
                header += " & ";
        }
        //latexString << R"(@{}} \hline)" << std::endl;
        latexString << R"(@{}})" << std::endl;
        latexString << header << R"(\\)" << std::endl;

        latexString << std::setprecision(m_config.yields_table_num_prec_yields);

        // Start with signals
        if (!signal_processes.empty()) {
            //latexString << "Signal sample" << ((signal_processes.size() == 1) ? "" : "s") << R"( & \\ \hline)" << std::endl;
            latexString << R"(\hline)" << std::endl;

            // Loop
            for (const auto& p: signal_processes_sorted) {

                latexString << p.second << " & & ";

                for (const auto& c: categories) {
                    std::string categ = c.second;

                    //latexString << "$" << signal_yields[categ][p].first << R"( {\scriptstyle\ \pm\ )" << std::sqrt(signal_yields[categ][p].second + std::pow(process_systematics[std::make_tuple(SIGNAL, categ, p)], 2)) << "}$ & ";
                    //latexString << "$" << signal_yields[categ][p].first << R"( {\scriptstyle\ \pm\ )" << std::sqrt(signal_yields[categ][p].second) << "}$ & ";
                    latexString << "$" << signal_yields[categ][p.first].first << "$ & $\\pm$ & $" << std::sqrt(signal_yields[categ][p.first].second) << "$ & ";
                }

                latexString.seekp(latexString.tellp() - 2l);
                latexString << R"( \\ )" << std::endl;
            }

            // Space
            //if (!mc_processes.empty() || has_data)
            //    latexString << R"( & \\)" << std::endl;
        }

        // Then MC samples
        if (!mc_processes.empty()) {
            //latexString << "SM sample" << ((mc_processes.size() == 1) ? "" : "s") << R"( & \\ \hline)" << std::endl;
            latexString << R"(\hline)" << std::endl;

            // Loop
            //for (const auto& p: mc_processes) {
            for (const auto& p: mc_processes_sorted) {

                latexString << p.second << " & & ";

                for (const auto& c: categories) {
                    std::string categ = c.second;

                    if (mc_yields[categ][p.first].first > 0)
                      //latexString << "$" << mc_yields[categ][p].first << R"( {\scriptstyle\ \pm\ )" << std::sqrt(mc_yields[categ][p].second + std::pow(process_systematics[std::make_tuple(MC, categ, p)], 2)) << "}$ & ";
                      //latexString << "$" << mc_yields[categ][p.first].first << R"( {\scriptstyle\ \pm\ )" << std::sqrt(mc_yields[categ][p.first].second) << "}$ & ";
                      latexString << "$" << mc_yields[categ][p.first].first << "$ &$\\pm$& $" << std::sqrt(mc_yields[categ][p.first].second) << "$ & ";
                    else
                      //latexString << "$" << "0" << R"( {\scriptstyle\ \pm\ )" << std::sqrt(mc_yields[categ][p.first].second) << "}$ & ";
                      latexString << "$" << "0" << "$ &$\\pm$& $" << std::sqrt(mc_yields[categ][p.first].second) << "$ & ";
                }

                latexString.seekp(latexString.tellp() - 2l);
                latexString << R"( \\ )" << std::endl;
            }

            // Space
            //latexString << R"( & \\ \hline)" << std::endl;
            //latexString << R"(Total {\scriptsize $\pm$ (stat) $\pm$ (syst)} & )";
            latexString << R"(Total {\scriptsize $\pm$ (stat)} & & )";

            for (const auto& c: categories) {
                //latexString << "$" << mc_total[c.second] << R"({\scriptstyle\ \pm\ )" << std::sqrt(mc_total_sqerrs[c.second]) << R"(\ \pm\ )" << std::sqrt(total_systematics_squared[c.second][MC]) << "}$ & ";
                //latexString << "$" << mc_total[c.second] << R"({\scriptstyle\ \pm\ )" << std::sqrt(mc_total_sqerrs[c.second]) << "}$ & ";
                latexString << "$" << mc_total[c.second] << "$ &$\\pm$& $" << std::sqrt(mc_total_sqerrs[c.second]) << "$ & ";
            }

            latexString.seekp(latexString.tellp() - 2l);
            latexString << R"( \\ )" << std::endl;
        }

        // Print data
        if (has_data) {
            latexString << R"(\hline)" << std::endl;
            //latexString << R"(Data {\scriptsize $\pm$ (stat)} & )";
            //latexString << R"(Data & )";
            latexString << R"(Data & & )";
            latexString << std::setprecision(0);

            for (const auto& c: categories) {
                // Compute poisson errors on the data yields
            //    static const double alpha = 1. - 0.682689492;
                int64_t yield = data_yields[c.second];
            //    double error_low = yield - ROOT::Math::gamma_quantile(alpha / 2., yield, 1.);
            //    double error_high = ROOT::Math::gamma_quantile_c(alpha / 2., yield, 1.) - yield;
            //    latexString << format_number_with_errors(yield, error_low, error_high, 0, m_config.yields_table_num_prec_yields) << " & ";
                latexString << yield << " & & & ";
            }

            latexString.seekp(latexString.tellp() - 2l);
            latexString << R"( \\ )" << std::endl;
        }

        // And finally data / MC
        if (!mc_processes.empty() && has_data) {
            latexString << R"(\hline)" << std::endl;
            latexString << R"(Data / Background prediction & & )";
            latexString << std::setprecision(m_config.yields_table_num_prec_ratio);

            for (const auto& c: categories) {
                std::string categ = c.second;
                int64_t data_yield = data_yields[categ];
                double ratio = data_yield / mc_total[categ];

                //static const double alpha = 1. - 0.682689492;
                //double error_data_low = data_yield - ROOT::Math::gamma_quantile(alpha / 2., data_yield, 1.);
                //double error_data_high = ROOT::Math::gamma_quantile_c(alpha / 2., data_yield, 1.) - data_yield;

                double error_mc = std::sqrt(mc_total_sqerrs[categ] + total_systematics_squared[categ][MC]);

                //double error_low = ratio * std::sqrt(std::pow(error_data_low / data_yields[categ], 2) +  std::pow(error_mc / mc_total[categ], 2));
                //double error_high = ratio * std::sqrt(std::pow(error_data_high / data_yields[categ], 2) +  std::pow(error_mc / mc_total[categ], 2));

                double error_low = std::pow(error_mc / mc_total[categ], 2);
                double error_high = std::pow(error_mc / mc_total[categ], 2);

                latexString << format_number_with_errors(ratio, error_low, error_high, m_config.yields_table_num_prec_ratio, m_config.yields_table_num_prec_ratio) << " & ";
                //latexString << ratio << " & ";
            }

            latexString.seekp(latexString.tellp() - 2l);
            latexString << R"( \\ )" << std::endl;
        }

        //latexString << R"(\hline)" << std::endl;
        latexString << R"(\end{tabular})" << std::endl;
    } else {
      std::cerr << "Error: yields table alignment " << m_config.yields_table_align << " is not recognized (for now, only \"h\" and \"v\" are supported)" << std::endl;
      return false;
    }

    if(CommandLineCfg::get().verbose)
      std::cout << "LaTeX yields table:\n\n" << latexString.str() << std::endl;

    fs::path outputName(m_outputPath);
    outputName /= "yields.tex";

    std::ofstream out(outputName.string());
    out << latexString.str();
    out.close();

    return true;
  }


  // systematics table
  bool plotIt::systematics(std::vector<Plot>::iterator plots_begin, std::vector<Plot>::iterator plots_end){
    std::cout << "Producing LaTeX systematic table.\n";

    std::map<std::string, double> data_yields;

    std::map< std::string, std::map<std::string, std::pair<double, double> > > mc_yields;
    std::map< std::string, double > mc_total;
    std::map< std::string, double > mc_total_sqerrs;
    std::set<std::string> mc_processes;

    std::map< std::string, std::map<std::string, std::pair<double, double> > > signal_yields;
    std::set<std::string> signal_processes;

    std::map<
        std::tuple<Type, std::string, std::string>, // Type, category, systematics name
        double
    > process_systematics;


    std::map<
        std::tuple<Type, std::string, std::string>, // Type, category, systematics name
        double
    > process_systematics_up;

    std::map<
        std::tuple<Type, std::string, std::string>, // Type, category, systematics name
        double
    > process_systematics_dn;

    std::map<
        std::string,
        std::map<Type, double>
    > total_systematics_squared;

    std::map<
        std::string,
        std::map<Type, double>
    > total_systematics_squared_up;

    std::map<
        std::string,
        std::map<Type, double>
    > total_systematics_squared_dn;

    std::vector< std::pair<int, std::string> > categories;

    bool has_data(false);

    for ( auto it = plots_begin; it != plots_end; ++it ) {
      auto& plot = *it;
      if (!plot.use_for_yields)
        continue;

      if (plot.yields_title.find("$") == std::string::npos)
          replace_substr(plot.yields_title, "_", "\\_");

      if( std::find_if(categories.begin(), categories.end(), [&](const std::pair<int, std::string> &x){ return x.second == plot.yields_title; }) != categories.end() )
          continue;
      categories.push_back( std::make_pair(plot.yields_table_order, plot.yields_title) );

      std::map<std::tuple<Type, std::string>, double> plot_total_systematics;
      std::map<std::tuple<Type, std::string>, double> plot_total_systematics_up;
      std::map<std::tuple<Type, std::string>, double> plot_total_systematics_dn;

      // Open all files, and find histogram in each
      for (auto& file: m_files) {
        if (! loadObject(file, plot)) {
          std::cout << "Could not retrieve plot from " << file.path << std::endl;
          return false;
        }

        if ( file.type == DATA ){
          TH1* h = dynamic_cast<TH1*>(file.object);
          data_yields[plot.yields_title] += h->Integral(0, h->GetNbinsX() + 1);
          has_data = true;
          continue;
        }

        std::string process_name = file.yields_group;

        if (process_name.find("$") == std::string::npos)
            replace_substr(process_name, "_", "\\_");

        if (process_name.find("#") != std::string::npos) {
            // We assume it's a ROOT LaTeX string. Enclose the string into $$, and replace
            // '#' by '\'

            replace_substr(process_name, "#", R"(\)");
            process_name = "$" + process_name + "$";
        }

        std::pair<double, double> yield_sqerror;
        TH1* hist( dynamic_cast<TH1*>(file.object) );

        if (m_config.generated_events_histogram.length() > 0 and file.generated_events < 0 and file.type != DATA ) {
            //generated_events = - 1.0 if not declared in file yaml
            std::shared_ptr<TFile> input(TFile::Open(file.path.c_str()));
            TH1* hevt = dynamic_cast<TH1*>(input->Get(m_config.generated_events_histogram.c_str()));
            file.generated_events = hevt->GetBinContent(m_config.generated_events_bin);
        }
        double factor = file.cross_section * file.branching_ratio / file.generated_events;

        if (! m_config.no_lumi_rescaling) {
          factor *= m_config.luminosity.at(file.era);
        }
        if (!CommandLineCfg::get().ignore_scales)
          factor *= m_config.scale * file.scale;

        if (!plot.is_rescaled)
          hist->Scale(factor);

        for (auto& syst: *file.systematics) {
          syst.update();
          syst.scale(factor);
        }

        // Retrieve yield and stat. error, taking overflow into account
        yield_sqerror.first = hist->IntegralAndError(0, hist->GetNbinsX() + 1, yield_sqerror.second);
        yield_sqerror.second = std::pow(yield_sqerror.second, 2);

        // Add systematics
        double file_total_systematics = 0;
        double file_total_systematics_up = 0;
        double file_total_systematics_dn = 0;
        for (auto& syst: *file.systematics) {

          TH1* nominal_shape = static_cast<TH1*>(syst.nominal_shape.get());
          TH1* up_shape = static_cast<TH1*>(syst.up_shape.get());
          TH1* down_shape = static_cast<TH1*>(syst.down_shape.get());

          if (! nominal_shape || ! up_shape || ! down_shape)
              continue;

          double nominal_integral = nominal_shape->Integral(0, nominal_shape->GetNbinsX() + 1);
          double up_integral = up_shape->Integral(0, up_shape->GetNbinsX() + 1);
          double down_integral = down_shape->Integral(0, down_shape->GetNbinsX() + 1);

          // For asym. error, define up/dn separately.
          // Be careful on the sign
          double temp_syst_error_up = up_integral - nominal_integral;
          double temp_syst_error_dn = down_integral - nominal_integral;

          double total_syst_error = std::max(
                  std::abs(up_integral - nominal_integral),
                  std::abs(nominal_integral - down_integral)
          );

          double total_syst_error_up;
          double total_syst_error_dn;

          // Normal case, double-sided
          if (temp_syst_error_up * temp_syst_error_dn <= 0) {
            if (temp_syst_error_up >= 0 and temp_syst_error_dn < 0) {
              total_syst_error_up = temp_syst_error_up;
              total_syst_error_dn = temp_syst_error_dn;
            } else {
              total_syst_error_up = temp_syst_error_dn;
              total_syst_error_dn = temp_syst_error_up;
            }
          }
          // One-sided.
          else {
            if (temp_syst_error_up > 0) {
                total_syst_error_up = std::max(temp_syst_error_up, temp_syst_error_dn);
                total_syst_error_dn = 0.0;
            } else {
                total_syst_error_up = 0.0;
                total_syst_error_dn = std::max(temp_syst_error_up, temp_syst_error_dn);
            }
          }

          file_total_systematics += total_syst_error * total_syst_error;
          file_total_systematics_up += total_syst_error_up * total_syst_error_up;
          file_total_systematics_dn += total_syst_error_dn * total_syst_error_dn;

          auto key = std::make_tuple(file.type, syst.name());
          plot_total_systematics[key] += total_syst_error;
          plot_total_systematics_up[key] += total_syst_error_up;
          plot_total_systematics_dn[key] += total_syst_error_dn;
        }

        // file_total_systematics contains the quadratic sum of all the systematics for this file
        process_systematics[std::make_tuple(file.type, plot.yields_title, process_name)] += std::sqrt(file_total_systematics);
        process_systematics_up[std::make_tuple(file.type, plot.yields_title, process_name)] += std::sqrt(file_total_systematics_up);
        process_systematics_dn[std::make_tuple(file.type, plot.yields_title, process_name)] += std::sqrt(file_total_systematics_dn);

        if ( file.type == MC ){
          ADD_PAIRS(mc_yields[plot.yields_title][process_name], yield_sqerror);
          mc_total[plot.yields_title] += yield_sqerror.first;
          mc_total_sqerrs[plot.yields_title] += yield_sqerror.second;
          mc_processes.emplace(process_name);
        }
        if ( file.type == SIGNAL ){
          ADD_PAIRS(signal_yields[plot.yields_title][process_name], yield_sqerror);
          signal_processes.emplace(process_name);
        }
      }

      // Get the total systematics for this category
      for (auto& syst: plot_total_systematics) {
        total_systematics_squared[plot.yields_title][std::get<0>(syst.first)] += syst.second * syst.second;
      }
      for (auto& syst: plot_total_systematics_up) {
        total_systematics_squared_up[plot.yields_title][std::get<0>(syst.first)] += syst.second * syst.second;
      }
      for (auto& syst: plot_total_systematics_dn) {
        total_systematics_squared_dn[plot.yields_title][std::get<0>(syst.first)] += syst.second * syst.second;
      }
    }

    if( ( !(mc_processes.size()+signal_processes.size()) && !has_data ) || !categories.size() ){
      std::cout << "No processes/data/categories defined\n";
      return false;
    }

    //std::vector<std::string> mc_processes_sorted;
    std::map< std::string, std::string > mc_processes_sorted;
    for(auto &process_name: mc_processes) {
      auto process_name_modified = process_name;
      if (process_name.find_first_of("0123456789") != std::string::npos) {
          // trick to sort process: put integer in front of the process name
          // then we remove here
          // Assume we have < 10 processes

          auto begin = process_name.find_first_of("0123456789");
          process_name_modified = process_name.substr(begin+1);

          std::string mathmark = "$";
          if (process_name.compare(0, 1, mathmark) == 0) process_name_modified = "$" + process_name_modified;
      }
      //mc_processes_sorted.emplace_back(process_name_modified);
      mc_processes_sorted[process_name] = process_name_modified;
    }

    std::map< std::string, std::string > signal_processes_sorted;
    for(auto &process_name: signal_processes) {
      auto process_name_modified = process_name;
      if (process_name.find_first_of("0123456789") != std::string::npos) {

          auto begin = process_name.find_first_of("0123456789");
          process_name_modified = process_name.substr(begin+1);

          std::string mathmark = "$";
          if (process_name.compare(0, 1, mathmark) == 0) process_name_modified = "$" + process_name_modified;
      }
      signal_processes_sorted[process_name] = process_name_modified;
    }

    // Sort according to user-defined order
    std::sort(categories.begin(), categories.end(), [](const std::pair<int, std::string>& cat1, const std::pair<int, std::string>& cat2){  return cat1.first < cat2.first; });

    std::ostringstream latexString;
    latexString << std::setiosflags(std::ios_base::fixed);
    if (m_config.yields_table_align.find("v") != std::string::npos) {

    // Start with signals
    if (!signal_processes.empty()) {
      for (const auto& p: signal_processes) {
        latexString << p << " & ";

        for (const auto& c: categories) {
          std::string categ = c.second;
          //latexString << R"($\pm$ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(std::pow(process_systematics[std::make_tuple(SIGNAL, categ, p)], 2))) << R"(\% & )";
          //latexString << R"($\pm$ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(std::pow(process_systematics[std::make_tuple(SIGNAL, categ, p)], 2))/signal_yields[categ][p].first)*100 << R"(\% & )";
          latexString << R"($^{+)"\
                      << std::setprecision(m_config.yields_table_num_prec_yields)\
                      << (std::sqrt(std::pow(process_systematics_up[std::make_tuple(SIGNAL, categ, p)], 2))/signal_yields[categ][p].first)*100\
                      <<  R"(\%}_{-)"\
                      << (std::sqrt(std::pow(process_systematics_dn[std::make_tuple(SIGNAL, categ, p)], 2))/signal_yields[categ][p].first)*100\
                      << R"(\%}$ & )";
        }

        latexString.seekp(latexString.tellp() - 2l);
        latexString << R"( \\ )" << std::endl;
      }
    }

    latexString << R"(\hline)" << std::endl;

    // Then MC samples
    if (!mc_processes.empty()) {
      for (const auto& p: mc_processes_sorted) {
          latexString << p.second << " & ";

          for (const auto& c: categories) {
            std::string categ = c.second;
            if(mc_yields[categ][p.first].first > 0) {
              //latexString << R"($\pm$ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(std::pow(process_systematics[std::make_tuple(MC, categ, p)], 2))) << R"(\% & )";
              //latexString << R"($\pm$ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(std::pow(process_systematics[std::make_tuple(MC, categ, p)], 2))/mc_yields[categ][p].first)*100 << R"(\% & )";
              latexString << R"($^{+)"\
                          << std::setprecision(m_config.yields_table_num_prec_yields)\
                          << (std::sqrt(std::pow(process_systematics_up[std::make_tuple(MC, categ, p.first)], 2))/mc_yields[categ][p.first].first)*100
                          <<  R"(\%}_{-)"\
                          << (std::sqrt(std::pow(process_systematics_dn[std::make_tuple(MC, categ, p.first)], 2))/mc_yields[categ][p.first].first)*100
                          << R"(\%}$ & )";
            } else {
              latexString << R"($\pm$ )" << R"($-$ )" << R"(\% & )";
            }
          }

          latexString.seekp(latexString.tellp() - 2l);
          latexString << R"( \\ )" << std::endl;
        }

        latexString << R"(\hline)" << std::endl;
        latexString << R"(Total sys. unc. & )";

        for (const auto& c: categories) {
          //latexString << R"($\pm $ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(total_systematics_squared[c.second][MC])) << R"(\% & )";
          //latexString << R"($\pm $ )" << std::setprecision(m_config.yields_table_num_prec_yields) << (std::sqrt(total_systematics_squared[c.second][MC])/mc_total[c.second])*100 << R"(\% & )";
          latexString << R"($^{+)"
                      << std::setprecision(m_config.yields_table_num_prec_yields)\
                      << (std::sqrt(total_systematics_squared_up[c.second][MC])/mc_total[c.second])*100\
                      <<  R"(\%}_{-)"\
                      << (std::sqrt(total_systematics_squared_dn[c.second][MC])/mc_total[c.second])*100\
                      << R"(\%}$ & )";
        }
        latexString.seekp(latexString.tellp() - 2l);
        latexString << R"( \\ )" << std::endl;
      }

    } else {
      std::cerr << "Error: systematics table alignment " << m_config.yields_table_align << " is not recognized (for now, only \"h\" and \"v\" are supported)" << std::endl;
      return false;
    }

    if(CommandLineCfg::get().verbose)
      std::cout << "LaTeX systematics table:\n\n" << latexString.str() << std::endl;

    fs::path outputName(m_outputPath);
    outputName /= "systematics.tex";

    std::ofstream out(outputName.string());
    out << latexString.str();
    out.close();

    return true;
  }

  void plotIt::plotAll() {

    m_style.reset(createStyle(m_config));

    // First, explode plots to match all glob patterns

    std::vector<Plot> plots;
    if (m_config.mode == "tree") {
      plots = m_plots;
    } else {
      if (!expandObjects(m_files[0], plots)) {
        return;
      }
    }

    if (!m_config.book_keeping_file_name.empty()) {
      fs::path outputName = m_outputPath / m_config.book_keeping_file_name;
      m_config.book_keeping_file.reset(TFile::Open(outputName.native().c_str(), "recreate"));
    }

    constexpr std::size_t plots_per_chunk = 20;

    auto plots_begin = plots.begin();
    auto plots_end = plots.begin();
    while ( plots_end != plots.end() ) {
      plots_begin = plots_end;
      if ( std::distance(plots_begin, plots.end()) > plots_per_chunk ) {
        plots_end = plots_begin+plots_per_chunk;
      } else {
        plots_end = plots.end();
      }

      if (CommandLineCfg::get().verbose)
          std::cout << "Loading plots " << std::distance(plots.begin(), plots_begin) << "-" << std::distance(plots.begin(), plots_end) << " of " << plots.size() << "..." << std::endl;

      for (File& file: m_files) {
        if (! loadAllObjects(file, plots_begin, plots_end))
            return;
      }

      if (CommandLineCfg::get().verbose)
          std::cout << "done." << std::endl;

      if (CommandLineCfg::get().do_plots) {
        for ( auto it = plots_begin; it != plots_end; ++it ) {
          plotIt::plot(*it);
        }
      }

      if (CommandLineCfg::get().do_yields) {
        plotIt::yields(plots_begin, plots_end);
      }

      if (CommandLineCfg::get().do_systematics) {
        plotIt::systematics(plots_begin, plots_end);
      }
    }

    for (File& file: m_files) {
      file.handle.reset();
      file.friend_handles.clear();
    }

    if (m_config.book_keeping_file) {
      m_config.book_keeping_file->Close();
      m_config.book_keeping_file.reset();
    }
  }

  bool plotIt::loadAllObjects(File& file, std::vector<Plot>::const_iterator plots_begin, std::vector<Plot>::const_iterator plots_end) {

    file.object = nullptr;
    file.objects.clear();

    if (m_config.mode == "tree") {

        if (!file.chain.get()) {
          file.chain.reset(new TChain(m_config.tree_name.c_str()));
          file.chain->Add(file.path.c_str());
        }

        for ( auto it = plots_begin; it != plots_end; ++it ) {
          const auto& plot = *it;

          auto x_axis_range = plot.log_x ? plot.log_x_axis_range : plot.x_axis_range;

          std::shared_ptr<TH1> hist(new TH1F((plot.uid + std::to_string(file.id)).c_str(), "", plot.binning_x, x_axis_range.start, x_axis_range.end));
          hist->SetDirectory(gROOT);

          file.chain->Draw((plot.draw_string + ">>" + plot.uid + std::to_string(file.id)).c_str(), plot.selection_string.c_str());

          hist->SetDirectory(nullptr);
          
          file.objects.emplace(plot.uid, hist.get());

          TemporaryPool::get().addRuntime(hist);
        }

        return true;
    }

    if (! file.handle)
      file.handle.reset(TFile::Open(file.path.c_str()));
    if (! file.handle)
      return false;

    file.systematics_cache.clear();

    for ( auto it = plots_begin; it != plots_end; ++it ) {
      const auto& plot = *it;

      std::string plot_name = plot.name;

      // Rename plot name according to user's transformations
      plot_name = applyRenaming(file.renaming_ops, plot_name);

      TObject* obj = file.handle->Get(plot_name.c_str());

      if (obj) {
        std::shared_ptr<TObject> cloned_obj(obj->Clone());

        // Rebin by array
        if (std::static_pointer_cast<TH1>(cloned_obj) && std::size(plot.rebin_arr) > 0) {
            auto h = std::static_pointer_cast<TH1>(cloned_obj).get();
            auto size_arr = std::size(plot.rebin_arr);
            double xbins[size_arr];
            for (size_t i=0; i<size_arr; i++) {
                xbins[i] = plot.rebin_arr.at(i);
            }
            auto hnew = h->Rebin(int(size_arr)-1, plot_name.c_str(), xbins);
            cloned_obj.reset(hnew);
        }


        TemporaryPool::get().addRuntime(cloned_obj);

        file.objects.emplace(plot.uid, cloned_obj.get());

        if (file.type != DATA) {
          for (auto& syst: m_systematics) {
              if (std::regex_search(file.path, syst->on))
                  file.systematics_cache[plot.uid].push_back(syst->newSet(cloned_obj.get(), file, plot));
          }
          for (auto& syst: m_systematics_siglike) {
              if (std::regex_search(file.path, syst->on))
                  file.systematics_cache_siglike[plot.uid].push_back(syst->newSet(cloned_obj.get(), file, plot));
          }
        }

        continue;
      }

      std::cout << "Error: object '" << plot_name << "' inheriting from '" << plot.inherits_from << "' not found in file '" << file.path << "'" << std::endl;
      return false;
    }

    return true;
  }

  bool plotIt::loadObject(File& file, const Plot& plot) {

    file.object = nullptr;

    auto it = file.objects.find(plot.uid);

    if (it == file.objects.end()) {
      auto exception = std::runtime_error("Object not found in cache. It should be here since it was preloaded before. Object name: " + plot.name + " in " + file.path);
      std::cerr << exception.what() << std::endl;
      throw exception;
    }

    file.object = it->second;

    file.systematics = & file.systematics_cache[plot.uid];
    file.systematics_siglike = & file.systematics_cache_siglike[plot.uid];

    return true;
  }

  bool plotIt::expandFiles() {
    std::vector<File> files;

    for (File& file: m_files) {
      std::vector<std::string> matchedFiles = glob(file.path);
      if (matchedFiles.empty()) {
          std::cerr << "Error: no files matching '" << file.path << "' (either the file does not exist, or the expression does not match any file)" << std::endl;
          return false;
      }
      for (std::string& matchedFile: matchedFiles) {
        File f = file;
        f.path = matchedFile;

        files.push_back(f);
      }
    }

    m_files = files;

    return true;
  }

  /**
   * Merge the labels of the global configuration and the current plot.
   * If some are duplicated, only keep the plot label
   **/
  std::vector<Label> plotIt::mergeLabels(const std::vector<Label>& plotLabels) {
    std::vector<Label> labels = plotLabels;

    // Add labels from global configuration, and check for duplicates
    for (auto& globalLabel: m_config.labels) {

      bool duplicated = false;
      for (auto& label: plotLabels) {
        if (globalLabel.text == label.text) {
          duplicated = true;
          break;
        }
      }

      if (! duplicated)
        labels.push_back(globalLabel);
    }

    return labels;
  }

  void get_directory_content(TDirectory* root, const std::string& prefix, std::vector<std::string>& content) {
      TIter it(root->GetListOfKeys());
      TKey* key = nullptr;

      while ((key = static_cast<TKey*>(it()))) {
          std::string name = key->GetName();
          std::string cl = key->GetClassName();

          if (cl.find("TDirectory") != std::string::npos) {
              std::string new_prefix = prefix;
              if (!prefix.empty())
                  new_prefix += "/";
              new_prefix += name;
              get_directory_content(static_cast<TDirectory*>(key->ReadObj()), new_prefix, content);
          } else if (cl.find("TH") != std::string::npos) {
              if (name.find("__") != std::string::npos) {
                  // TODO: Maybe we should be a bit less strict and check that the
                  // systematics specified is included in the configuration file?
                  continue;
              }

              std::string new_prefix = prefix;
              if (!new_prefix.empty())
                  new_prefix += "/";
              content.push_back(new_prefix + name);
          }
      }
  }

  /**
   * Open 'file', and expand all plots
   */
  bool plotIt::expandObjects(File& file, std::vector<Plot>& plots) {
    file.object = nullptr;
    plots.clear();

    // Optimization. Look if any of the plots have a glob pattern (either *, ? or [)
    // If not, do not iterate of the file to match pattern, it's useless
    std::vector<Plot> glob_plots;
    for (Plot& plot: m_plots) {
        if ((plot.name.find("*") != std::string::npos) || (plot.name.find("?") != std::string::npos) || (plot.name.find("[") != std::string::npos)) {
            glob_plots.push_back(plot);
        } else {
            plots.push_back(plot.Clone(plot.name));
        }
    }

    if (glob_plots.empty()) {
        return true;
    }

    std::shared_ptr<TFile> input(TFile::Open(file.path.c_str()));
    if (! input.get())
      return false;

    // Create file structure, flattening any directory
    TIter root_keys(input->GetListOfKeys());

    std::vector<std::string> file_content;
    get_directory_content(input.get(), "", file_content);

    for (Plot& plot: glob_plots) {
        bool match = false;
        std::vector<std::string> matched;

        for (const auto& content: file_content) {

            // Check name
            if (fnmatch(plot.name.c_str(), content.c_str(), FNM_CASEFOLD) == 0) {

                // Check if this name is excluded
                if ((plot.exclude.length() > 0) && (fnmatch(plot.exclude.c_str(), content.c_str(), FNM_CASEFOLD) == 0)) {
                    continue;
                }

                // The same object can be stored multiple time with a different key
                // The iterator returns first the object with the highest key, which is the most recent object
                // Check if we already have a plot with the same exact name
                if (std::find_if(matched.begin(), matched.end(), [&content](const std::string& p) { return p == content; }) != matched.end()) {
                    continue;
                }

                // Got it!
                match = true;
                matched.push_back(content);
                plots.push_back(plot.Clone(content));
            }
        }

        if (! match) {
            std::cout << "Warning: object '" << plot.name << "' inheriting from '" << plot.inherits_from << "' does not match something in file '" << file.path << "'" << std::endl;
        }
    }

    if (!plots.size()) {
      std::cout << "Error: no plots found in file '" << file.path << "'" << std::endl;
      return false;
    }

    return true;
  }

  std::shared_ptr<PlotStyle> plotIt::getPlotStyle(const File& file) {
    if (file.legend_group.length() && m_legend_groups.count(file.legend_group)) {
      return m_legend_groups[file.legend_group].plot_style;
    } else {
      return file.plot_style;
    }
  }
}

int main(int argc, char** argv) {

  try {

    TCLAP::CmdLine cmd("Plot histograms", ' ', "0.1");

    TCLAP::ValueArg<std::string> histogramsFolderArg("i", "histograms-folder", "histograms base folder (default: current directory)", false, "./", "string", cmd);

    TCLAP::ValueArg<std::string> outputFolderArg("o", "output-folder", "output folder", true, "", "string", cmd);

    TCLAP::ValueArg<std::string> eraArg("e", "era", "era to restrict to", false, "", "string", cmd);

    TCLAP::SwitchArg ignoreScaleArg("", "ignore-scales", "Ignore any scales present in the configuration file", cmd, false);

    TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose output (print summary)", cmd, false);

    TCLAP::SwitchArg yieldsArg("y", "yields", "Produce LaTeX table of yields", cmd, false);

    TCLAP::SwitchArg systematicsArg("s", "systematics", "Produce LaTeX table of systematics", cmd , false);

    TCLAP::SwitchArg plotsArg("p", "plots", "Do not produce the plots - can be useful if only the yields table is needed", cmd, false);

    TCLAP::SwitchArg unblindArg("u", "unblind", "Unblind the plots, ie ignore any blinded-range in the configuration", cmd, false);

    TCLAP::SwitchArg systematicsBreakdownArg("b", "systs-breadown", "Print systematics details for each MC process separately in addition to the total contribution", cmd, false);

    TCLAP::UnlabeledValueArg<std::string> configFileArg("configFile", "configuration file", true, "", "string", cmd);

    TCLAP::SwitchArg qcdArg("q", "qcd", "Process with QCD samples (histo name with QCD)", cmd, false);

    TCLAP::SwitchArg dyArg("d", "dyincl", "Process with DY inclusive sample (DYJetsToLL_M50_amc)", cmd, false);

    TCLAP::SwitchArg allSigArg("a", "allSig", "plot only vector signal for LFV, default for ", cmd, false);

    TCLAP::SwitchArg noSigArg("n", "noSig", "plot no signal for LFV, default for ", cmd, false);

    TCLAP::ValueArg<std::string> selectSigArg("", "selectSig", "plot only select signal for LFV", false, "", "string", cmd);

    TCLAP::SwitchArg desytopArg("t", "desytop", "Switch for DESY top framework", cmd, false);

    TCLAP::SwitchArg binnedArg("", "binned", "Draw with 'ttbarsignal' binned samples", cmd, false);

    cmd.parse(argc, argv);

    //bool isData = dataArg.isSet();

    fs::path histogramsPath(fs::canonical(histogramsFolderArg.getValue()));

    if (! fs::exists(histogramsPath)) {
      std::cout << "Error: histograms path " << histogramsPath << " does not exist" << std::endl;
    }

    fs::path outputPath(outputFolderArg.getValue());

    if (! fs::exists(outputPath)) {
      std::cout << "Error: output path " << outputPath << " does not exist" << std::endl;
      return 1;
    }

    if( plotsArg.getValue() && !yieldsArg.getValue() ) {
      std::cerr << "Error: we have nothing to do" << std::endl;
      return 1;
    }

    CommandLineCfg::get().era = eraArg.getValue();
    CommandLineCfg::get().ignore_scales = ignoreScaleArg.getValue();
    CommandLineCfg::get().verbose = verboseArg.getValue();
    CommandLineCfg::get().do_plots = !plotsArg.getValue();
    CommandLineCfg::get().do_yields = yieldsArg.getValue();
    CommandLineCfg::get().do_systematics = systematicsArg.getValue();
    CommandLineCfg::get().unblind = unblindArg.getValue();
    CommandLineCfg::get().systematicsBreakdown = systematicsBreakdownArg.getValue();
    CommandLineCfg::get().do_qcd = qcdArg.getValue();
    CommandLineCfg::get().dyincl = dyArg.getValue();
    CommandLineCfg::get().allSig = allSigArg.getValue();
    CommandLineCfg::get().noSig = noSigArg.getValue();
    CommandLineCfg::get().selectSig = selectSigArg.getValue();
    CommandLineCfg::get().desytop = desytopArg.getValue();
    CommandLineCfg::get().binned = binnedArg.getValue();

    plotIt::plotIt p(outputPath);
    if (!p.parseConfigurationFile(configFileArg.getValue(), histogramsPath))
        return 1;

    p.plotAll();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
