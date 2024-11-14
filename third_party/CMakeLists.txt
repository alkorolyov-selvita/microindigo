if (BUILD_STANDALONE)
    # InChI don't have a Conan package yet
    add_subdirectory(cppcodec)
    add_subdirectory(inchi)
    add_subdirectory(libpng)
    add_subdirectory(rapidjson)
    add_subdirectory(zlib)
    add_subdirectory(tinyxml2)
endif()
if (ENABLE_TESTS)
    add_subdirectory(googletest)
endif()

add_subdirectory(object_threadsafe)
