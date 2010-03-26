#include "itkImplicitFunctionBase.h"
#include "itkObjectFactory.h"

namespace itk
{
template< typename TOutput, unsigned int VPointDimension, class TInput >
class ImplicitFunctionBaseTestHelper : 
  public ImplicitFunctionBase< TOutput, VPointDimension, TInput >
{
  public:
    typedef ImplicitFunctionBaseTestHelper  Self;
    typedef ImplicitFunctionBase< TOutput, VPointDimension, TInput > 
                                            Superclass;
    typedef SmartPointer< Self >            Pointer;
    typedef SmartPointer< const Self >      ConstPointer;

    itkStaticConstMacro( PointDimension, unsigned int,
        VPointDimension );

    /** \brief Input type for the ImplicitFunctionBase. */
    typedef typename Superclass::InputType InputType;

    /** \brief Output type for the ImplicitFunctionBase value. */
    typedef typename Superclass::OutputType OutputType;

    /** \brief Define type vector for gradient vector. */
    typedef typename Superclass::GradientType GradientType;

    /** \brief Define type vector for hessian matrix. */
    typedef typename Superclass::HessianType HessianType;

    itkNewMacro(Self);

    itkTypeMacro(ImplicitFunctionBaseTestHelper,ImplicitFunctionBase);

    OutputType Evaluate( const InputType& iPt )
      { return static_cast< OutputType >( 0 ); }

    GradientType Gradient( const InputType& iPt )
      { 
      GradientType oG; 
      oG.Fill( static_cast< OutputType >( 0 ) );
      return oG;
      }

    HessianType Hessian( const InputType& iPt )
      {
      HessianType oH;
      oH.Fill( static_cast< OutputType >( 0 ) );
      return oH;
      }

  protected:
    ImplicitFunctionBaseTestHelper() {}
    ~ImplicitFunctionBaseTestHelper() {}
  private:
    ImplicitFunctionBaseTestHelper( const Self& );
    void operator = ( const Self& );
};
}


int main( int argc, char* argv[] )
{
  typedef itk::Point< float, 2 > PointType;
  typedef itk::ImplicitFunctionBaseTestHelper< double, 2, PointType > 
    FunctionType;
  FunctionType::Pointer f = FunctionType::New();

  PointType p;
  p.Fill( 0. );

  f->Evaluate( p );
  f->Gradient( p );
  f->Hessian( p );
  f->Laplacian( p );

  return EXIT_SUCCESS;
}
