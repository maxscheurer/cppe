#pragma once
#include <sstream>
#include <string>
#include <vector>

template <typename Out>
void split(const std::string& s, char delim, Out result);

std::vector<std::string> split(const std::string& s, char delim);

std::string trim(const std::string& str, const std::string& whitespace = " \t");

std::string reduce(const std::string& str, const std::string& fill = " ",
                   const std::string& whitespace = " \t\n");
