PROJECT( itkGenLS )

CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

FIND_PACKAGE( ITK )

IF( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ENDIF( ITK_FOUND )

INCLUDE_DIRECTORIES( BEFORE 
  ${itkGenLS_SOURCE_DIR}/src/
#  ${itkGenLS_SOURCE_DIR}/utilities/
  ${itkGenLS_BINARY_DIR}/src/
#  ${itkGenLS_BINARY_DIR}/utilities/
)

INCLUDE( CTest )

SET( EXECUTABLE_OUTPUT_PATH 
   ${itkGenLS_BINARY_DIR}/bin 
   CACHE PATH "Single output directory for building all executables." 
 ) 
  
SET( LIBRARY_OUTPUT_PATH 
  ${itkGenLS_BINARY_DIR}/lib 
  CACHE PATH "Single output directory for building all libraries." 
)


#ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( test )


