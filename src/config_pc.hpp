//
// Created by 69029 on 6/25/2021.
//

#ifndef UNLAYERED_LIBRA_CONFIG_PC_HPP
#define UNLAYERED_LIBRA_CONFIG_PC_HPP

#ifdef USE_VIRGO
#include <virgo/polyCommit.hpp>
#define F   virgo::fieldElement
#define F_ONE   virgo::fieldElement::one()
#define F_ZERO  virgo::fieldElement::zero()
#endif

#ifdef USE_HYRAX_P224
#include <hyrax-p224/src/polyVerifier.hpp>
#define F   hyrax_p224::fieldElement
#define G   hyrax_p224::groupElement
#define F_ONE   hyrax_p224::fieldElement::one()
#define F_ZERO  hyrax_p224::fieldElement::zero()
#endif

#endif //UNLAYERED_LIBRA_CONFIG_PC_HPP
