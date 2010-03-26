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
      virtual OutputType Evaluate( const InputType& iPt ) = 0;

      /**
      *  \brief Evaluate the gradient of the ImplicitFunctionBase at a given position.
      *  iPt is transformed through mTransform (if provided).
      *
      *  \param iPt
      *  \return
      */
      virtual GradientType Gradient( const InputType& iPt ) = 0;


      GradientType Normal( const InputType& iPt );

      /**
      *  \brief Evaluate the hessian of the ImplicitFunctionBase at a given position.
      *
      *  \param iPt
      *  \return
      */
      virtual HessianType Hessian( const InputType& iPt ) = 0;

      /**
      * \brief returns the Laplacian Value at iPt.
      *
      * \param[in] iPt
      * \return \f[\sum_{i=0}^{dimension}\frac{\partial^2 f}{\partial x_i^2}\f]
      */
      OutputType Laplacian( const InputType& iPt );

      /**
      *  \brief returns in the same time the value and the gradient of the
      *  implicit ImplicitFunctionBase at a given Position.
      *
      *  \param[in] iPt
      */
      virtual void ValueAndGradient( const InputType& iPt,
         OutputType& oValue, GradientType& oGrad );

      /**
      *  \brief returns in the same time value, gradient and hessian of the
      *  implicit ImplicitFunctionBase at a given iPt.
      *
      *  \param[in] iPt
      *  \param[in] oValue
      *  \param[in] oGrad
      *  \param[in] oHess
      */
      virtual void ValueGradientAndHessian( const InputType& iPt,
         OutputType& oValue, GradientType& oGrad, HessianType& oHess );

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
      virtual OutputType TaubinErrorValue(
        const InputType& iPt,
        const OutputType& iThreshold = vnl_math::eps );
  

      bool IsInside( const InputType& iPt );

      bool IsOutside( const InputType& iPt );

      OutputType MeanCurvature( const InputType& iPt );

      OutputType GaussKroneckerCurvature( const InputType& iPt );

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
         GradientType& oEmax );

      HessianType CurvatureTensor( const InputType& iPt );


   protected:
      /** Default constructor */
      ImplicitFunctionBase( );

      /** Destructor */
      virtual ~ImplicitFunctionBase( );

      /** */
      virtual void PrintSelf(std::ostream& os, Indent indent) const;

      struct DispatchBase {};

      template< unsigned int >
      struct Dispatch : public DispatchBase {};

      OutputType MeanCurvature( const Dispatch< 3 >&,
        const InputType& iPt );

      OutputType MeanCurvature( const DispatchBase&,
        const InputType& iPt );

      OutputType GaussKroneckerCurvature( const Dispatch< 2 >&,
        const InputType& iPt );

      OutputType GaussKroneckerCurvature( const DispatchBase&,
        const InputType& iPt );

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
        const InputType& iPt );

      HessianType CurvatureTensor( const DispatchBase&,
        const InputType& iPt );

      void PrincipalCurvatures( const DispatchBase&,
        const InputType& iPt,
        OutputType& oKmin,
        OutputType& oKmax,
        GradientType& oEmin,
        GradientType& oEmax );

      void PrincipalCurvatures( const Dispatch<3>&,
        const InputType& iPt,
        OutputType& oKmin,
        OutputType& oKmax,
        GradientType& oEmin,
        GradientType& oEmax );

  private:
      /** Not implemented */
      ImplicitFunctionBase( const Self& );

      /** Not implemented */
      void operator=( const Self& );

   };

} // End of namespace itk

#include "itkImplicitFunctionBase.txx"

#endif
