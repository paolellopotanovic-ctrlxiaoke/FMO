/*
 * integral_operator_type.hpp
 *
 * Created: RELEASE_DATE
 * Author: Bryce M. Westheimer
 *
 */

#ifndef LIBACCINT_INCLUDE_INTEGRAL_INTEGRAL_OPERATOR_TYPE_HPP_
#define LIBACCINT_INCLUDE_INTEGRAL_INTEGRAL_OPERATOR_TYPE_HPP_

// External header(s)

// Internal header(s)

// STL header(s)

namespace libaccint {

	namespace host {

		namespace engine {

			/*! \enum IntegralOperatorType
			*   \ingoup TODO:ADD_GROUP
			*
			*   \brief Integral operator type enumerator for supported integrals
			*
			*   \details TODO
			*
			*   \author Bryce M. Westheimer
			*
			*   \date RELEASE_DATE
			*
			*/
			enum class IntegralOperatorType {
				Overlap,
				KineticEnergy,
				ElectronNuclearCoulomb,
				CoreHamiltonian,
				ElectronElectronCoulomb
				// TODO: Add as more are implemented
			};

		} // namespace libaccint::host::engine

	} // namespace libaccint::host

} // namespace libaccint

#endif /* LIBACCINT_INCLUDE_INTEGRAL_INTEGRAL_OPERATOR_TYPE_HPP_ */