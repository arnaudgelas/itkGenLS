#ifndef __itkImplicitPlaneFunction_h
#define __itkImplicitPlaneFunction_h

#include "itkObjectFactory.h"
#include "itkImplicitFunctionBase.h"
#include "itkFixedArray.h"

namespace itk
{
  template< typename TOutput,
    unsigned int VPointDimension,
    class TInput >
  class ImplicitPlaneFunction :
    public ImplicitFunctionBase< TOutput, VPointDimension, TInput >
  {
  public:
    typedef ImplicitPlaneFunction Self;
    typedef ImplicitFunctionBase< TOutput, VPointDimension, TInput > 
      Superclass;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    itkNewMacro( Self );

    itkTypeMacro( ImplicitPlaneFunction, ImplicitFunctionBase );

    itkStaticConstMacro( PointDimension, unsigned int,
      VPointDimension );

    typedef typename Superclass::InputType InputType;
    typedef typename Superclass::OutputType OutputType;
    typedef typename Superclass::GradientType GradientType;
    typedef typename Superclass::HessianType HessianType;

    typedef FixedArray< OutputType, PointDimension+1 > CoefficientVectorType;

    itkSetMacro( Coefficients, CoefficientVectorType );
    itkGetMacro( Coefficients, CoefficientVectorType );

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

#include "itkImplicitPlaneFunction.txx"
#endif
