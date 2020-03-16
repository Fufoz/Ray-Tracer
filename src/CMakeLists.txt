
add_executable(hoth 
	main.cc
	logging.cc
	bitmap.cc
)

target_link_libraries(hoth fmt)
target_include_directories(hoth PRIVATE ${CMAKE_SOURCE_DIR}/extern/stb)
install(TARGETS hoth RUNTIME DESTINATION ./)