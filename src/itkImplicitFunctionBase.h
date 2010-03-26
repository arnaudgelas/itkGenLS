/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-10

 Copyright (c) 2009-10, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef __itkImplicitFunctionBase_h
#define __itkImplicitFunctionBase_h

#include "itkObject.h"
#include "itkVector.h"
#include "itkMatrix.h"

#include "vnl/vnl_math.h"
// #include <vnl/vnl_matrix.h>

 // Ron Goldman "curvature formulas for implicit curves and surfaces"
 // Torsten Langer "Curvature Computations for Implicit Submanifolds"
namespace itk
{
   /**
   *  \defgroup implicit ImplicitFunctionBase
   */

   /**
   *  \class ImplicitFunctionBase
   *
   *  \brief \ingroup implicit
   *  Abstract class which represents an implicit ImplicitFunctionBase.
   *  Any point <em>p(x,y,z)</em> has an associate scalar value <em>f(p)</em>.
   *
   *  \note According to the litterature, if
   *        \li <em> f(p) > 0   p</em>  is outside,
   *        \li <em> f(p) = 0   p</em> lays in the surface,
   *        \li <em> f(p) < 0   p</em> is inside.
   *
   */

   template < typename TOutput,
      unsigned int VPointDimension,
      class TInput >
  class ITK_EXPORT ImplicitFunctionBase : public Object
  {
  public:

      /** Standard class typedefs.      */
      typedef ImplicitFunctionBase Self;

      typedef Object Superclass;

      typedef SmartPointer< Self >  Pointer ;
      typedef SmartPointer< const Self >  ConstPointer ;

      /** Run-time type information (and related methods). */
      itkTypeMacro( ImplicitFunctionBase, Object ) ;

      itkStaticConstMacro( PointDimension, unsigned int,
        VPointDimension );

      /** \brief Input type for the ImplicitFunctionBase. */
      typedef TInput InputType;

      /** \brief Output type for the ImplicitFunctionBase value. */
      typedef TOutput OutputType;

      /** \brief Define type vector for gradient vector. */
      typedef Vector< OutputType, PointDimension > GradientType;

      /** \brief Define type vector for hessian matrix. */
      typedef Matrix< OutputType, PointDimension, PointDimension > HessianType;

      /**
      *  \brief Evaluate the ImplicitFunctionBase at a given position.
      *
      *  \param
      *  \return
      */
      virtual OutputType Evaluate( const InputType& iPt ) const = 0;

      /**
      *  \brief Evaluate the gradient of the ImplicitFunctionBase at a given position.
      *  iPt is transformed through mTransform (if provided).
      *
      *  \param iPt
      *  \return
      */
      virtual GradientType Gradient( const InputType& iPt ) const = 0;
      GradientType Normal( const InputType& iPt ) const
        {
        GradientType grad = Gradient( iPt );
        OutputType sq_norm = grad.GetSquaredNorm();
        if( sq_norm > vnl_math::eps )
          {
          return grad / vcl_sqrt( sq_norm );
          }
        else
          {
          return grad;
          }
        }

      /**
      *  \brief Evaluate the hessian of the ImplicitFunctionBase at a given position.
      *
      *  \param iPt
      *  \return
      */
      virtual HessianType Hessian( const InputType& iPt ) const = 0;

      /**
      * \brief returns the Laplacian Value at iPt.
      *
      * \param[in] iPt
      * \return \f[\sum_{i=0}^{dimension}\frac{\partial^2 f}{\partial x_i^2}\f]
      */
      OutputType Laplacian( const InputType& iPt ) const
        {
        HessianType hess = Hessian( iPt );

        OutputType oValue = 0.;

        for( unsigned int i = 0; i < PointDimension; ++i )
          {
          oValue += hess[i][i];
          }
        return oValue;
        }

      /**
      *  \brief returns in the same time the value and the gradient of the
      *  implicit ImplicitFunctionBase at a given Position.
      *
      *  \param[in] iPt
      */
      virtual void ValueAndGradient( const InputType& iPt,
         OutputType& oValue, GradientType& oGrad ) const
        {
          oValue = Evaluate( iPt );
          oGrad = Gradient( iPt );
        }

      /**
      *  \brief returns in the same time value, gradient and hessian of the
      *  implicit ImplicitFunctionBase at a given iPt.
      *
      *  \param[in] iPt
      *  \param[in] oValue
      *  \param[in] oGrad
      *  \param[in] oHess
      */
      void ValueGradientAndHessian( const InputType& iPt,
         OutputType& oValue, GradientType& oGrad, HessianType& oHess ) const
        {
          ValueAndGradient( iPt, oValue, oGrad );
          oHess = Hessian( iPt );
        }

      /**
      *  \brief Evaluate Approximation of the distance (\e Dist)
      *  between one point and the 0-isosurface.
      *
      *  \f[ Dist = \frac{\left|f(p)\right|}
      *     {\left\|\vec{\nabla}f(p)\right\|}\f]
      *
      *  @param iPt
      *  @return = Dist
      *     if(\f$ \left\|\vec{\nabla}f(p)\right\|\neq 0\f$)
      *  @return = \f$ \left|f(p)\right| \f$
      *     else
      */
      virtual OutputType ImplicitFunctionBaseTaubinErrorValue( const InputType& iPt,
        const OutputType& iThreshold = vnl_math::eps ) const
        {
        OutputType oValue = Evaluate( iPt );

        GradientType grad = Gradient( iPt );

        if( oValue != 0. )
          {
          // Implicit ImplicitFunctionBase value
          oValue = vnl_math_abs( oValue );

          // Norm of Gradient Vector
          OutputType GradNorm = grad.GetNorm( );

          if( GradNorm < iThreshold )
            {
            GradNorm = iThreshold;
            }

          if( GradNorm > vnl_math::eps )
            {
            oValue /= GradNorm;
            }
          }

        return oValue;
        }

      bool IsInside( const InputType& iPt ) const
        {
        OutputType value = Evaluate( iPt );
        return ( value < 0. );
        }

      bool IsOutside( const InputType& iPt ) const
        {
        OutputType value = Evaluate( iPt );
        return ( value > 0. );
        }

      OutputType MeanCurvature( const InputType& iPt ) const
        {
        return MeanCurvature( Dispatch< PointDimension >(), iPt );
        }

      OutputType GaussKroneckerCurvature( const InputType& iPt ) const
        {
        return GaussKroneckerCurvature( Dispatch< PointDimension >(), iPt );
        }

      /**
      *  \brief Compute the Curvatures Tensor
      *
      *  @param iPt
      *  @param oKmin \f$\kappa_{min}\f$
      *  @param oKmax \f$\kappa_{max}\f$
      *  @param oEmin corresponding principal direction for \f$\kappa_{min}\f$.
      *  @param oEmax corresponding principal direction for \f$\kappa_{max}\f$.
      */
      void PrincipalCurvaturesTensor( const InputType& iPt,
         OutputType& oKmin,
         OutputType& oKmax,
         GradientType& oEmin,
         GradientType& oEmax ) const
        {

        PrincipalCurvaturesTensor( Dispatch< PointDimension >(), iPt,
          oKmin, oKmax, oEmin, oEmax );
        }

      HessianType CurvatureTensor( const InputType& iPt ) const
        {
        return CurvatureTensor( Dispatch< PointDimension >(), iPt );
        }


   protected:
      /** Default constructor */
      ImplicitFunctionBase( )  : Superclass() {}

      /** Destructor */
      virtual ~ImplicitFunctionBase( ) {}

      /** */
      virtual void PrintSelf(std::ostream& os, Indent indent) const {};

      struct DispatchBase {};

      template< unsigned int >
      struct Dispatch : public DispatchBase {};

      OutputType MeanCurvature( const Dispatch< 3 >&,
        const InputType& iPt ) const
        {
        GradientType grad = Gradient( iPt );
        HessianType hessian = Hessian( iPt );

        OutputType norm_g = grad.GetNorm();

        if( norm_g == 0. )
          {
          return 0.;
          }
        else
          {
          HessianType Identity;
          Identity.SetIdentity();
          OutputType inv_norm_g = 1. / norm_g;
          grad *= inv_norm_g;
          OutputType factor = ( Identity - grad * grad.transpose() ) * hessian;

          return -inv_norm_g * factor;
          }
        }

      OutputType MeanCurvature( const DispatchBase&,
        const InputType& iPt ) const
        {
        GradientType grad = Gradient( iPt );
        HessianType  hessian = Hessian( iPt );

        OutputType sq_norm_g = grad.GetSquaredNorm();
        OutputType norm_g = vcl_sqrt( sq_norm_g );
        OutputType gHgt = grad * hessian * grad.transpose();
        OutputType trace_H = 0.;

        for( unsigned int dim = 0; dim < PointDimension; ++dim )
          {
          trace_H += hessian[dim][dim];
          }

        return ( 1. / ( ( PointDimension - 1 ) * norm_g * sq_norm_g ) *
          (gHgt - sq_norm_g * trace_H) );
        }

      OutputType GaussKroneckerCurvature( const Dispatch< 2 >&,
        const InputType& iPt ) const
        {
        return 0.;
        }

      OutputType GaussKroneckerCurvature( const DispatchBase&,
        const InputType& iPt ) const
        {
        GradientType g = Gradient( iPt );
        OutputType sq_norm_g = g.GetSquaredNorm();
        if( sq_norm_g == 0. )
          {
          return 0.;
          }
        else
          {
          OutputType norm_g = vcl_sqrt( sq_norm_g );
          HessianType H = Hessian( iPt );
          OutputType coeff = ( PointDimension % 2 == 1 ) ? -norm_g : norm_g;

          for( unsigned int dim = 0; dim < PointDimension; ++dim )
            {
            coeff *= norm_g;
            }
            // should be conjuguate but H is real...
          OutputType gtH_g = g.transpose() * H.transpose() * g;

          return gtH_g / coeff;
          }
        }

      /**
      *  \brief Compute the curvature tensor for a given position.
      *
      *  @param iPt
      *
      *  @return oTensor
      *  \f[ \frac{1}{\|\vec{\nabla}f\|}\biggl(\mathcal{H}_{f}-
      *     \frac{\vec{\nabla}f^{t}.\mathcal{H}_{f}.\vec{\nabla}f}
      *     {\|\vec{\nabla}f\|^{2}} \biggr)\f]
      *
      *
      *  \todo check the formula!
      */
      HessianType CurvatureTensor( const Dispatch<3>&,
        const InputType& iPt ) const
        {
        GradientType g = Gradient( iPt );
        HessianType H = Hessian( iPt );

        OutputType inv_sq_norm = 1. / g.GetSquaredNorm();
        OutputType inv_norm = vcl_sqrt( inv_sq_norm );

        Matrix< OutputType, PointDimension, 1 > Gradient;

        for( unsigned int i = 0; i < PointDimension; i++ )
            Gradient[i] = g[i];

        ImplicitFunctionBaseGradient2( iPt, Hessian );

        HessianType oTensor = ( H -
          inv_sq_norm * g.transpose() * H * g ) * inv_norm;

        return oTensor;
        }

      HessianType CurvatureTensor( const DispatchBase&,
        const InputType& iPt ) const
        {
        HessianType oTensor;
        oTensor.Fill( 0. );
        return oTensor;
        }

      void PrincipalCurvatures( const DispatchBase&,
        const InputType& iPt,
        OutputType& oKmin,
        OutputType& oKmax,
        GradientType& oEmin,
        GradientType& oEmax )
        {
        oKmin = 0.;
        oKmax = 0.;
        oEmin.Fill( 0. );
        oEmax.Fill( 0. );
        }

      void PrincipalCurvatures( const Dispatch<3>&,
        const InputType& iPt,
        OutputType& oKmin,
        OutputType& oKmax,
        GradientType& oEmin,
        GradientType& oEmax )
        {
        OutputType mean_curvature = MeanCurvature( Dispatch<3>(), iPt );
        OutputType gauss_curvature =
          GaussKroneckerCurvature( Dispatch<3>(), iPt );

        OutputType delta =
          vnl_math_abs( mean_curvature * mean_curvature - gauss_curvature );

        OutputType sqrt_delta = vcl_sqrt( delta );

        oKmax = mean_curvature + sqrt_delta;
        oKmin = mean_curvature - sqrt_delta;

        if( sqrt_delta < 1e-6 )
          {
          oEmin.Fill( 0. );
          oEmax.Fill( 0. );
          }
        else
          {
          HessianType tensor = CurvatureTensor( Dispatch<3>(), iPt );

          /* solve for eigenvectors */
          OutputType t[3];
          t[0] = tensor[0][0] + oKmax;
          t[1] = tensor[1][1] + oKmax;
          t[2] = tensor[2][2] + oKmax;

          GradientType v[3];
          OutputType len[3];

          v[0][0] = tensor[0][1] * tensor[1][2] - tensor[0][2] * t[1];
          v[1][0] = tensor[0][2] * tensor[1][0] - tensor[1][2] * t[0];
          v[2][0] = t[0] * t[1] - tensor[0][1] * tensor[1][0];
          len[0] = vcl_sqrt( v[0][0] * v[0][0] +
            v[1][0] * v[1][0] + v[2][0] * v[2][0] );

          v[0][1] = tensor[0][1] * t[2] - tensor[0][2] * tensor[2][1];
          v[1][1] = tensor[0][2] * tensor[2][0]- t[0] * t[2];
          v[2][1] = t[0] * tensor[2][1] - tensor[0][1] * tensor[2][0];
          len[1] = vcl_sqrt( v[0][1] * v[0][1] +
            v[1][1] * v[1][1] + v[2][1] * v[2][1] );

          v[0][2] = t[1] * t[2] - tensor[1][2] * tensor[2][1];
          v[1][2] = tensor[1][2] * tensor[2][0]- tensor[1][0] * t[2];
          v[2][2] = tensor[1][0] * tensor[2][1] - tensor[2][0]* t[1];
          len[2] = vcl_sqrt( v[0][2] * v[0][2] +
            v[1][2] * v[1][2] + v[2][2] * v[2][2] );

          unsigned int index = 0;
          OutputType len_max = len[index];

          if( len[1] > len_max )
            {
            index = 1;
            len_max = len[1];
            }
          if( len[2] > len_max )
            {
            index = 2;
            len_max = len[2];
            }

          oEmax = v[index] / len[index];

          GradientType n = Normal( iPt );

          oEmin[0] = oEmax[1] * n[2] - oEmax[2] * n[1];
          oEmin[1] = oEmax[2] * n[0] - oEmax[0] * n[2];
          oEmin[2] = oEmax[0] * n[1] - oEmax[1] * n[0];
          }
        }

  private:
      /** Not implemented */
      ImplicitFunctionBase( const Self& );

      /** Not implemented */
      void operator=( const Self& );

   };

} // End of namespace itk

#endif
