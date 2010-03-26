SET( TEST_SRC 
  ImplicitFunctionBaseTest
  )

FOREACH( var ${TEST_SRC} )
  ADD_EXECUTABLE( ${var} ${var}.cxx )
  TARGET_LINK_LIBRARIES( ${var} ITKCommon )
ENDFOREACH( var ${TEST_SRC} )

ADD_TEST( ImplicitFunctionBaseTest
  ${EXECUTABLE_OUTPUT_PATH}/ImplicitFunctionBaseTest 
)
