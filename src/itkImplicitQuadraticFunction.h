#ifndef __itkImplicitQuadraticFunction_h
#define __itkImplicitQuadraticFunction_h

#include "itkObjectFactory.h"
#include "itkImplicitFunctionBase.h"
#include "itkFixedArray.h"

namespace itk
{
  template< typename TOutput,
    unsigned int VPointDimension,
    class TInput >
  class ImplicitQuadraticFunction :
    public ImplicitFunctionBase< TOutput, VPointDimension, TInput >
  {
  public:
    typedef ImplicitQuadraticFunction Self;
    typedef ImplicitFunctionBase< TOutput, VPointDimension, TInput > 
      Superclass;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    itkNewMacro( Self );

    itkTypeMacro( ImplicitQuadraticFunction, ImplicitFunctionBase );

    itkStaticConstMacro( PointDimension, unsigned int,
      VPointDimension );

    typedef typename Superclass::InputType InputType;
    typedef typename Superclass::OutputType OutputType;
    typedef typename Superclass::GradientType GradientType;
    typedef typename Superclass::HessianType HessianType;

    typedef FixedArray< OutputType, 
      PointDimension * ( PointDimension + 1 ) / 2 +
      PointDimension + 1 
    > CoefficientVectorType;

    itkSetMacro( Coefficients, CoefficientVectorType );
    itkGetMacro( Coefficients, CoefficientVectorType );

    OutputType Evaluate( const InputType& iPt );
    GradientType Gradient( const InputType& iPt );
    HessianType Hessian( const InputType& iPt );

  protected:
    ImplicitQuadraticFunction();
    ~ImplicitQuadraticFunction();

    CoefficientVectorType m_Coefficients;

  private:
    ImplicitQuadraticFunction( const Self& );
    void operator = ( const Self& );
  };
}

#include "itkImplicitQuadraticFunction.txx"
#endif
