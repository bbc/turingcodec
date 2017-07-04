cmake_minimum_required(VERSION 2.8)

find_package(Git)
if(GIT_FOUND)
    include(${CMAKE_CURRENT_SOURCE_DIR}/GetGitRevisionDescription.cmake)
    git_describe(GIT_DESCRIBE)
else()
    set(GIT_DESCRIBE version-unknown)
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/git-describe.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/git-describe.h)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/..")

if (USE_SYSTEM_BOOST)
else()
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../boost")
endif()

if(MSVC)
    #We statically link to reduce dependancies
    foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        endif(${flag_var} MATCHES "/MD")
        if(${flag_var} MATCHES "/MDd")
            string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
        endif(${flag_var} MATCHES "/MDd")
    endforeach(flag_var)
endif(MSVC)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_ALL_NO_LIB /bigobj /GS- /MP")
endif()

set(LINK_LIBRARIES boost_program_options boost_filesystem boost_system boost_timer boost_chrono)

if (CMAKE_COMPILER_IS_GNUCC)
    add_compile_options(-fno-operator-names -pthread -std=c++1y)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-strict-aliasing")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -DVALGRIND_FRIENDLY")
    set(LINK_LIBRARIES ${LINK_LIBRARIES} -pthread rt)
endif (CMAKE_COMPILER_IS_GNUCC)

set(SEI_HEADERS
    sei/ReadSei.h
    sei/active_parameter_sets.h
    sei/all.h
    sei/alpha_channel_info.h
    sei/bsp_initial_arrival_time.h
    sei/bsp_nesting.h
    sei/buffering_period.h
    sei/chroma_resampling_filter_hint.h
    sei/colour_remapping_info.h
    sei/content_light_level.h
    sei/decoded_picture_hash.h
    sei/decoding_unit_info.h
    sei/deinterlaced_field_identification.h
    sei/display_orientation.h
    sei/filler_payload.h
    sei/film_grain_characteristics.h
    sei/frame_field_info.h
    sei/frame_packing_arrangement.h
    sei/inter_layer_constrained_tile_sets.h
    sei/knee_function_info.h
    sei/layers_not_present.h
    sei/mastering_display_colour_volume.h
    sei/no_display.h
    sei/overlay_info.h
    sei/pan_scan_rect.h
    sei/pic_timing.h
    sei/picture_snapshot.h
    sei/post_filter_hint.h
    sei/progressive_refinement_segment_end.h
    sei/progressive_refinement_segment_start.h
    sei/recovery_point.h
    sei/region_refresh_info.h
    sei/reserved_sei_message.h
    sei/scalable_nesting.h
    sei/scene_info.h
    sei/segmented_rect_frame_packing_arrangement.h
    sei/structure_of_pictures_info.h
    sei/temporal_motion_constrained_tile_sets.h
    sei/temporal_mv_prediction_constraints.h
    sei/temporal_sub_layer_zero_index.h
    sei/time_code.h
    sei/tone_mapping_info.h
    sei/user_data_registered_itu_t_t35.h
    sei/user_data_unregistered.h)

add_library(turing
    Cabac.cpp
    ContextModel.cpp
    EncodeResidual.cpp
    EncodeResidualPopcnt.cpp
    EncodeResidual.hpp
    EncodeResidual.h
    encode.cpp
    psnr.cpp
    TaskEncodeInput.cpp
    decode.cpp
    Encoder.cpp
    Rdoq.cpp
    TaskEncodeOutput.cpp
    DecodeCtu.cpp
    InputQueue.cpp
    Reconstruct.cpp
    TaskEncodeSubstream.cpp
    DecodeRbsp.cpp
    ThreadPool.cpp
    md5.c
    ScalingMatrices.cpp
    Picture.cpp
    Search.cpp
    SearchLzcnt.cpp
    sao.cpp
    signature.cpp
    ProgressReporter.cpp
    ScanOrder.cpp
    TaskDeblock.cpp
    testdecode.cpp
    TaskSao.cpp
    Aps.h
    Binarization.h
    BitField.h
    BlockData.h
    Cabac.h
    CabacWriter.h
    CandModeList.h
    CodeCache.h
    CodedData.h
    ConformanceStreams.h
    ContextModel.h
    Cost.h
    Decode.h
    DecodedPicture.h
    Dsp.h
    Encoder.h
    EncSao.h
    EstimateRate.h
    EstimateIntraComplexity.h
    FixedPoint.h
    GlobalState.h
    Handlers.h
    HevcMath.h
    HevcTypes.h
    InputQueue.h
    IntraReferenceSamples.h
    Levels.h
    LoopFilter.h
    md5.h
    Measure.h
    Memory.h
    MotionVector.h
    Mvp.h
    Padding.h
    Read.h
    Read.hpp
    Picture.h
    Global.h
    Primitive.h
    Profiler.h
    Profiles.h
    ProgressReporter.h
    RangeLimits.h
    RateControl.h
    RateControl.cpp
    Rdoq.h
    Reconstruct.h
    ReconstructionCache.h
    sao.h
    ScalingMatrices.h
    ScanOrder.h
    SCDetection.h
    Search.h
    Search.hpp
    Snake.h
    Speed.h
    StateCollocatedMotion.h
    StateDecode.h
    StateEncode.h
    StateFunctionTables.h
    StateParameterSets.h
    StatePicture.h
    StatePictures.h
    StateSpatial.h
    StateValues.h
    StateWavefront.h
    StreamReader.h
    StructMacros.h
    Syntax.h
    SyntaxCtu.hpp
    SyntaxElements.h
    SyntaxNal.hpp
    SyntaxRbsp.hpp
    SyntaxSei.h
    TaskDeblock.h
    TaskSao.h
    TaskEncodeInput.h
    TaskEncodeOutput.h
    TaskEncodeSubstream.h
    ThreadPool.h
    turing.h
    Violation.h
    Write.h
    ${SEI_HEADERS})

if (CMAKE_COMPILER_IS_GNUCC)
    set_source_files_properties(EncodeResidualPopcnt.cpp PROPERTIES COMPILE_FLAGS "-mlzcnt -mpopcnt")
    set_source_files_properties(SearchLzcnt.cpp PROPERTIES COMPILE_FLAGS "-mlzcnt -mpopcnt")
endif (CMAKE_COMPILER_IS_GNUCC)

source_group(SEI FILES ${SEI_HEADERS})

target_link_libraries (turing LINK_PUBLIC havoc)

unset(CMAKE_IMPORT_LIBRARY_SUFFIX)

add_executable (turing-exe main.cpp)

set_target_properties(turing-exe PROPERTIES OUTPUT_NAME turing)

target_link_libraries (turing-exe LINK_PUBLIC turing ${LINK_LIBRARIES})

if(UNIX)
	include(FindPkgConfig QUIET)
	if(PKG_CONFIG_FOUND)
		configure_file("libturing.pc.in" "libturing.pc" @ONLY)
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libturing.pc" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../boost/libs/chrono/src/libboost_chrono.a" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/boost")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../boost/libs/filesystem/src/libboost_filesystem.a" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/boost")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../boost/libs/program_options/src/libboost_program_options.a" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/boost")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../boost/libs/system/src/libboost_system.a" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/boost")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../boost/libs/timer/src/libboost_timer.a" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/boost")
	endif()
endif()

install(TARGETS turing LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES turing.h DESTINATION include)
