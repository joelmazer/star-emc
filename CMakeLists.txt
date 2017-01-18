cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(star-base)


# Add to path in order to pick up the FindXXX.cmake files included in this project
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/star-cmake")

# This project depends on ROOT
find_package(ROOT)


if(NOT ROOT_FOUND)
	message(FATAL_ERROR "Fatal error: ROOT package not found")
endif()

# Perform some setup standard to STAR experiment environment
include(StarCommon)

add_definitions(-D__ROOT__)


include_directories(
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${STAR_BASE_DIR}
	${STAR_STEVENT_DIR}
	${STAR_STDB_DIR}
	${STAR_DAQ_DIR}
	${STAR_MUDST_DIR}
	${STAR_PAMS_DIR}/global/inc
	${STAR_TABLES_DIR}
	${STAR_TABLES_DIR}/tables
)


# StEmcUtil library
root_generate_linkdef_and_dictionary(StEmcUtil HEADERS
	StEmcUtil/others/StEmcMath.h
	StEmcUtil/filters/StEmcFilter.h
	StEmcUtil/voltageCalib/VoltCalibrator.h
	StEmcUtil/voltageCalib/GainVoltCoeffCalculator.h
	StEmcUtil/database/StBemcTables.h
	StEmcUtil/database/StEmcMappingDb.h
	StEmcUtil/database/StEmcDecoder.h
	StEmcUtil/database/StEmcDbHandler.h
	StEmcUtil/database/StBemcTablesWriter.h
	StEmcUtil/geometry/StEmcGeom.h
	StEmcUtil/projection/StEmcPosition.h
)

add_library(StEmcUtil SHARED
	StEmcUtil/others/StEmcMath.cxx
	StEmcUtil/filters/StEmcFilter.cxx
	StEmcUtil/voltageCalib/VoltCalibrator.cxx
	StEmcUtil/voltageCalib/PmtIdentifier.cxx
	StEmcUtil/voltageCalib/GainVoltPmtParameters.cxx
	StEmcUtil/voltageCalib/GainVoltCoeffCalculator.cxx
	StEmcUtil/database/StEmcDecoder.cxx
	StEmcUtil/database/StBemcTables.cxx
	StEmcUtil/database/StEmcMappingDb.cxx
	StEmcUtil/database/StBemcTablesWriter.cxx
	StEmcUtil/database/StEmcDbHandler.cxx
	StEmcUtil/geometry/StEmcGeom.cxx
	StEmcUtil/projection/StEmcPosition.cxx
	StEmcUtil_dict.cxx
)


# StEEmcUtil library
root_generate_linkdef_and_dictionary(StEEmcUtil HEADERS
	StEEmcUtil/EEevent/EEeventDst.h
	StEEmcUtil/EEevent/EEsectorDst.h
	StEEmcUtil/EEevent/EEsmdHitDst.h
	StEEmcUtil/EEevent/EEtwHitDst.h
	StEEmcUtil/database/EEmcDbCrate.h
	StEEmcUtil/database/EEmcDbItem.h
	StEEmcUtil/database/StEEmcDb.h
	StEEmcUtil/EEmcGeom/EEmcGeomSimple.h
	StEEmcUtil/EEfeeRaw/EEfeeRawEvent.h
	StEEmcUtil/EEfeeRaw/EEmcEventHeader.h
	StEEmcUtil/EEfeeRaw/EEfeeDataBlock.h
	StEEmcUtil/EEfeeRaw/EEmcL3Tracks.h
	StEEmcUtil/EEfeeRaw/EEmcL3Tracks.h
	StEEmcUtil/EEfeeRaw/BsmdRawData.h
	StEEmcUtil/EEfeeRaw/EEstarTrig.h
)

add_library(StEEmcUtil SHARED
	StEEmcUtil/EEevent/EEeventDst.cxx
	StEEmcUtil/EEevent/EEsectorDst.cxx
	StEEmcUtil/EEevent/EEsmdHitDst.cxx
	StEEmcUtil/EEevent/EEtwHitDst.cxx
	StEEmcUtil/database/EEmcDbCrate.cxx
	StEEmcUtil/database/EEmcDbItem.cxx
	StEEmcUtil/database/StEEmcDb.cxx
	StEEmcUtil/EEmcGeom/EEmcGeomSimple.cxx
	StEEmcUtil/EEfeeRaw/BsmdRawData.cxx
	StEEmcUtil/EEfeeRaw/EEfeeDataBlock.cxx
	StEEmcUtil/EEfeeRaw/EEfeeRawEvent.cxx
	StEEmcUtil/EEfeeRaw/EEmcEventHeader.cxx
	StEEmcUtil/EEfeeRaw/EEmcL3Tracks.cxx
	StEEmcUtil/EEfeeRaw/EEname2Index.cxx
	StEEmcUtil/EEfeeRaw/EEstarTrig.cxx
	StEEmcUtil/EEfeeRaw/RootWrapper.cxx
	StEEmcUtil_dict.cxx
)


# A collective target to build all libraries in this project. Can be used to
# build all specified targets from a parent project
add_custom_target(star-emc DEPENDS StEmcUtil StEEmcUtil)


# Installation section
install(TARGETS StEmcUtil StEEmcUtil
	DESTINATION "${STAR_ADDITIONAL_INSTALL_PREFIX}/lib" OPTIONAL)
