#pragma once
#include <string>
#include <vector>

namespace libcppe {

struct version {
  /** Return the major part of the version */
  static int major_part();

  /** Return the minor part of the version */
  static int minor_part();

  /** Return the patch part of the version */
  static int patch_part();

  /** Is the compiled version a Debug version */
  static bool is_debug();

  /**  Return the version as a string */
  static std::string version_string();
};

/** Return the authors string */
std::string __authors__();

/** Return the contributors string */
std::string __contributors__();

/** Return the email string */
std::string __email__();

}  // namespace libcppe
