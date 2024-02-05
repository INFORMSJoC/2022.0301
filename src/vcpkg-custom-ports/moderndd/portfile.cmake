vcpkg_from_gitlab(
  GITLAB_URL
  https://gitlab.kuleuven.be
  OUT_SOURCE_PATH
  SOURCE_PATH
  REF
  v0.0.6
  SHA512
  3c4f86579f816acb67d346010ebff3afdbd0826641142024c6e84d45f12db27bed68836d4cd50b7fe2f7ef40618ec9e80f681899fb44c250db0466d58be0e331
  REPO
  u0056096/DecisionDiagramsModernCpp
)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME "ModernDD")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")
vcpkg_copy_pdbs()
