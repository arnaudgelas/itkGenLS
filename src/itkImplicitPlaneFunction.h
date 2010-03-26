#ifndef __itkImplicitPlaneFunction_h
#define __itkImplicitPlaneFunction_h

#include "itkImplicitFunctionBase.h"

namespace itk
{
  template< typename TOutput,
    unsigned int VPointDimension,
    class TInput >
  class ImplicitPlaneFunction :
    public ImplicitFunctionBase< TOutput, VPointDimension, TInput >
  {
  public:
    typedef FixedArray< OutputType, PointDimension+1 > CoefficientVectorType;

    OutputType Evaluate( const InputType& iPt );
    GradientType Gradient( const InputType& iPt );
    HessianType Hessian( const InputType& iPt );

  protected:
    ImplicitPlaneFunction();
    ~ImplicitPlaneFunction();

    CoefficientVectorType m_Coefficients;

  private:
    ImplicitPlaneFunction( const Self& );
    void operator = ( const Self& );
  };
}
#endif
