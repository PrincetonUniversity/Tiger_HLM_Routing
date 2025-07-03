#include "config_loader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

// Implementation of SimpleYamlParser methods

void SimpleYamlParser::parseFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();
    
    parseLines(lines);
}

std::string SimpleYamlParser::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

std::string SimpleYamlParser::removeQuotes(const std::string& str) {
    std::string trimmed = trim(str);
    
    // Remove inline comments first
    size_t commentPos = trimmed.find('#');
    if (commentPos != std::string::npos) {
        trimmed = trim(trimmed.substr(0, commentPos));
    }
    
    // Then remove quotes
    if (!trimmed.empty() && 
        ((trimmed.front() == '"' && trimmed.back() == '"') ||
         (trimmed.front() == '\'' && trimmed.back() == '\''))) {
        return trimmed.substr(1, trimmed.length() - 2);
    }
    return trimmed;
}

bool SimpleYamlParser::isInlineArray(const std::string& str) {
    std::string trimmed = trim(str);
    return !trimmed.empty() && trimmed.front() == '[' && trimmed.back() == ']';
}

std::vector<std::string> SimpleYamlParser::parseInlineArray(const std::string& str) {
    std::vector<std::string> result;
    std::string trimmed = trim(str);
    
    if (!isInlineArray(trimmed)) {
        return result;
    }
    
    // Remove brackets
    std::string content = trimmed.substr(1, trimmed.length() - 2);
    
    // Split by comma
    std::stringstream ss(content);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        std::string cleanItem = trim(item);
        if (!cleanItem.empty()) {
            result.push_back(removeQuotes(cleanItem));
        }
    }
    
    return result;
}

int SimpleYamlParser::getIndentLevel(const std::string& line) {
    int indent = 0;
    for (char c : line) {
        if (c == ' ') indent++;
        else if (c == '\t') indent += 4; // Treat tab as 4 spaces
        else break;
    }
    return indent;
}

bool SimpleYamlParser::isArrayItem(const std::string& line) {
    std::string trimmed = trim(line);
    return !trimmed.empty() && trimmed[0] == '-';
}

bool SimpleYamlParser::isComment(const std::string& line) {
    std::string trimmed = trim(line);
    return trimmed.empty() || trimmed[0] == '#';
}

void SimpleYamlParser::parseLines(const std::vector<std::string>& lines) {
    std::vector<std::string> sectionPath;
    std::string currentArrayKey;
    std::vector<std::map<std::string, std::string>> currentArray;
    std::vector<std::string> currentSimpleArray;
    bool inArray = false;
    bool inSimpleArray = false;
    int baseIndent = -1;
    
    for (size_t i = 0; i < lines.size(); i++) {
        const std::string& line = lines[i];
        
        if (isComment(line)) continue;
        
        int indent = getIndentLevel(line);
        std::string trimmed = trim(line);
        
        if (trimmed.empty()) continue;
        
        // Handle array items
        if (isArrayItem(trimmed)) {
            std::string arrayContent = trim(trimmed.substr(1)); // Remove '-'
            
            // Check if this is a simple array (just values) or complex array (objects)
            if (arrayContent.find(':') != std::string::npos) {
                // Complex array item (object)
                if (!inArray) {
                    inArray = true;
                    inSimpleArray = false;
                    currentArray.clear();
                }
                
                std::map<std::string, std::string> arrayItem;
                
                // Parse the current line
                size_t colonPos = arrayContent.find(':');
                if (colonPos != std::string::npos) {
                    std::string key = trim(arrayContent.substr(0, colonPos));
                    std::string value = removeQuotes(trim(arrayContent.substr(colonPos + 1)));
                    arrayItem[key] = value;
                }
                
                // Look ahead for more properties of this array item
                for (size_t j = i + 1; j < lines.size(); j++) {
                    std::string nextLine = lines[j];
                    if (isComment(nextLine)) continue;
                    
                    int nextIndent = getIndentLevel(nextLine);
                    std::string nextTrimmed = trim(nextLine);
                    
                    if (nextTrimmed.empty()) continue;
                    
                    // If next line is less indented or another array item, stop
                    if (nextIndent <= indent || isArrayItem(nextTrimmed)) {
                        break;
                    }
                    
                    // Parse key-value pair
                    size_t nextColonPos = nextTrimmed.find(':');
                    if (nextColonPos != std::string::npos) {
                        std::string nextKey = trim(nextTrimmed.substr(0, nextColonPos));
                        std::string nextValue = removeQuotes(trim(nextTrimmed.substr(nextColonPos + 1)));
                        arrayItem[nextKey] = nextValue;
                        i = j; // Skip these lines in main loop
                    }
                }
                
                currentArray.push_back(arrayItem);
            } else {
                // Simple array item (just a value)
                if (!inSimpleArray) {
                    inSimpleArray = true;
                    inArray = false;
                    currentSimpleArray.clear();
                }
                currentSimpleArray.push_back(removeQuotes(arrayContent));
            }
            continue;
        }
        
        // If we were in an array and now we're not, save it
        if ((inArray || inSimpleArray) && !isArrayItem(trimmed)) {
            if (inArray && !currentArray.empty()) {
                arrayMap[currentArrayKey] = currentArray;
            } else if (inSimpleArray && !currentSimpleArray.empty()) {
                simpleArrayMap[currentArrayKey] = currentSimpleArray;
            }
            inArray = false;
            inSimpleArray = false;
            currentArray.clear();
            currentSimpleArray.clear();
        }
        
        // Handle regular key-value pairs
        size_t colonPos = trimmed.find(':');
        if (colonPos != std::string::npos) {
            std::string key = trim(trimmed.substr(0, colonPos));
            std::string value = trim(trimmed.substr(colonPos + 1));
            
            // Remove inline comments from value before checking if it's empty
            std::string cleanValue = value;
            size_t commentPos = cleanValue.find('#');
            if (commentPos != std::string::npos) {
                cleanValue = trim(cleanValue.substr(0, commentPos));
            }
            
            // Adjust section path based on indentation
            if (baseIndent == -1) baseIndent = indent;
            
            int level = (indent - baseIndent) / 2; // Assuming 2-space indentation
            if (level < 0) level = 0;
            
            // Adjust section path
            if (level < sectionPath.size()) {
                sectionPath.resize(level);
            }
            
            if (cleanValue.empty()) {
                // This is a section header
                sectionPath.push_back(key);
                currentArrayKey = getSectionKey(sectionPath);
            } else {
                // This is a key-value pair
                std::vector<std::string> fullPath = sectionPath;
                fullPath.push_back(key);
                std::string fullKey = getSectionKey(fullPath);
                
                // Check if the value is an inline array
                if (isInlineArray(value)) {
                    std::vector<std::string> arrayItems = parseInlineArray(value);
                    simpleArrayMap[fullKey] = arrayItems;
                } else {
                    keyValueMap[fullKey] = removeQuotes(value);
                }
            }
        }
    }
    
    // Save any remaining array
    if (inArray && !currentArray.empty()) {
        arrayMap[currentArrayKey] = currentArray;
    } else if (inSimpleArray && !currentSimpleArray.empty()) {
        simpleArrayMap[currentArrayKey] = currentSimpleArray;
    }
}

std::string SimpleYamlParser::getSectionKey(const std::vector<std::string>& path) {
    std::string result;
    for (size_t i = 0; i < path.size(); i++) {
        if (i > 0) result += ".";
        result += path[i];
    }
    return result;
}

std::string SimpleYamlParser::getString(const std::string& key, const std::string& defaultValue) {
    auto it = keyValueMap.find(key);
    return (it != keyValueMap.end()) ? it->second : defaultValue;
}

int SimpleYamlParser::getInt(const std::string& key, int defaultValue) {
    auto it = keyValueMap.find(key);
    if (it != keyValueMap.end()) {
        return std::stoi(it->second);
    }
    return defaultValue;
}

double SimpleYamlParser::getDouble(const std::string& key, double defaultValue) {
    auto it = keyValueMap.find(key);
    if (it != keyValueMap.end()) {
        return std::stod(it->second);
    }
    return defaultValue;
}


// Implementation of ConfigLoader methods
ModelConfig ConfigLoader::loadConfig(const std::string& filename) {
    SimpleYamlParser parser;
    parser.parseFile(filename);
    
    ModelConfig config;
    
   
    // Load time period
    config.start_date = parser.getString("time.start_date");
    config.calendar = parser.getString("time.calendar");

    //Load solver/simulation parameters
    config.dt = parser.getDouble("solver.dt");
    config.rtol = parser.getDouble("solver.rtol");
    config.atol = parser.getDouble("solver.atol");
    config.simulation_resolution = parser.getInt("solver.simulation_resolution");
    
    // Load parameters
    config.parameters_file = parser.getString("parameters.filename");

    // Load initial conditions
    config.initial_conditions_flag = parser.getInt("initial_conditions.flag");
    config.initial_value = parser.getDouble("initial_conditions.initial_value");
    config.initial_conditions_filename = parser.getString("initial_conditions.filename");
    config.initial_conditions_varname = parser.getString("initial_conditions.varname");
    config.initial_conditions_id_varname = parser.getString("initial_conditions.id_varname");

    // Load boundary conditions
    config.boundary_conditions_flag = parser.getInt("boundary_conditions.flag");
    config.boundary_conditions_resolution = parser.getInt("boundary_conditions.resolution");
    config.boundary_conditions_filename = parser.getString("boundary_conditions.filename");     
    config.boundary_conditions_varname = parser.getString("boundary_conditions.varname");
    config.boundary_conditions_id_varname = parser.getString("boundary_conditions.id_varname"); 

    // Load reservoir routing
    config.reservoir_routing_flag = parser.getInt("reservoir.flag");
    config.reservoir_file = parser.getString("reservoir.res_file");
    
    // Load runoff parameters
    config.input_flag = parser.getInt("runoff.flag");
    config.runoff_resolution = parser.getInt("runoff.resolution");
    config.chunk_size = parser.getInt("runoff.chunk_size");
    config.runoff_path = parser.getString("runoff.folder");
    config.runoff_varname = parser.getString("runoff.varname");
    config.runoff_id_varname = parser.getString("runoff.id_varname");

    // Load output options
    config.output_flag = parser.getInt("output.flag");
    config.min_level = parser.getInt("output.level");    
    if(config.min_level < 1){
        std::cerr << "Warning: Minimum level for output is set to " << config.min_level 
                  << ". It should be at least 1. Setting it to 1." << std::endl;
        config.min_level = 1; // Ensure min_level is at least 1
    }
    config.link_list_filename = parser.getString("output.link_list_filename");
    config.series_filepath = parser.getString("output.series_filepath");
    config.snapshot_filepath = parser.getString("output.snapshot_filepath");

    return config;
}