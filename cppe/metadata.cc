#include "metadata.hh"
#include <sstream>
#include <vector>

namespace libcppe {
namespace {
static const std::string static_version_string = "0.2.1";

static const std::vector<std::string> version_split = [](const std::string& in) {
  std::vector<std::string> parts;
  std::stringstream ss(in);
  std::string item;
  while (std::getline(ss, item, '.')) parts.push_back(item);
  return parts;
}(static_version_string);

static int get_version_part(size_t part) {
  int ret;
  std::stringstream ss(version_split[part]);
  ss >> ret;
  return ret;
}
}  // namespace

int version::major_part() { return get_version_part(0); }
int version::minor_part() { return get_version_part(1); }
int version::patch_part() { return get_version_part(2); }
bool version::is_debug() {
#ifdef NDEBUG
  return false;
#else
  return true;
#endif  // NDEBUG
}

std::string version::version_string() { return static_version_string; }

std::string __authors__() { return "Maximilian Scheurer"; }

std::string __contributors__() { return "Peter Reinholdt, Michael F. Herbst, Lori A. Burns"; }

std::string __email__() { return "maximilian.scheurer@iwr.uni-heidelberg.de"; }

}  // namespace libcppe