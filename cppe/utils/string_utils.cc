#include "string_utils.hh"

template <typename Out>
void split(const std::string& s, char delim, Out result) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

std::string trim(const std::string& str, const std::string& whitespace) {
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos) return "";  // no content

  const auto strEnd   = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str, const std::string& fill,
                   const std::string& whitespace) {
  // trim first
  auto result = trim(str, whitespace);

  // replace sub ranges
  auto beginSpace = result.find_first_of(whitespace);
  while (beginSpace != std::string::npos) {
    const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
    const auto range    = endSpace - beginSpace;

    result.replace(beginSpace, range, fill);

    const auto newStart = beginSpace + fill.length();
    beginSpace          = result.find_first_of(whitespace, newStart);
  }

  return result;
}