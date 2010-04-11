#include "itkPoint.h"
#include "itkImplicitQuadraticFunction.h"

int main( int argc, char** argv )
{
  const unsigned int Dimension = 2;
  typedef itk::Point< double, Dimension > PointType;
  typedef float OutputType;

  typedef itk::ImplicitQuadraticFunction< OutputType, Dimension, PointType >
    QuadraticType;
  QuadraticType::Pointer f = QuadraticType::New();
 
  QuadraticType::CoefficientVectorType coeff;
  coeff[0] = 1.; // x*x
  coeff[1] = 1.; // x*y
  coeff[2] = 0.; // y*y
  coeff[3] = 0.; // x
  coeff[4] = 0.; // y
  coeff[5] = -1.;// 1
  
  f->SetCoefficients( coeff );
   
  PointType p;
  p[0] = 1.;
  p[1] = 0.; 

  if( f->Evaluate( p ) != 3 )
    {
    std::cerr <<"Evaluate(p) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  QuadraticType::GradientType g = f->Gradient( p );
  if( ( g[0] != 1 ) || ( g[1] != 1 ) )
    {
    std::cerr <<"Gradient( p ) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  QuadraticType::HessianType h = f->Hessian( p );
  if( ( h[0][0] != 0 ) || ( h[0][1] != 0 ) ||
      ( h[1][0] != 0 ) || ( h[1][1] != 0 ) )
    {
    std::cerr <<"Hessian( p ) Failure" <<std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

