add_executable(mhm driver.F90)
target_link_libraries(mhm PRIVATE mhm_lib)
install(TARGETS mhm RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
