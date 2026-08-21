/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Fatemeh Nouri
 * Copyright (c) 2026 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/form-factors/parametric-bhkmnr2026.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricBHKMNR2026Test :
    public TestCase
{
    public:
        ParametricBHKMNR2026Test() :
            TestCase("parametric_BHKMNR2026_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            {
                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                        =  0.13957039;
                p["0->pipi::s_0@BHKMNR2026"]           =  0;
                p["0->pipi::s_in@BHKMNR2026"]          =  0.85051;
                p["0->pipi::a_(+,1)^4@BHKMNR2026"]     = -0.56083194;
                p["0->pipi::a_(+,1)^5@BHKMNR2026"]     = -0.18059761;
                p["0->pipi::a_(+,1)^6@BHKMNR2026"]     = -0.31252685;
                p["0->pipi::a_(+,1)^7@BHKMNR2026"]     = -0.03797998;
                p["0->pipi::a_(+,1)^8@BHKMNR2026"]     = -0.1769615;
                p["0->pipi::a_(+,1)^9@BHKMNR2026"]     =  0;
                p["0->pipi::a_(+,1)^10@BHKMNR2026"]    =  0;
                p["0->pipi::a_(+,1)^11@BHKMNR2026"]    =  0;
                p["0->pipi::a_(+,1)^12@BHKMNR2026"]    =  0;
                p["0->pipi::M_(+,1,0)@BHKMNR2026"]     =  0.76033012;
                p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"] =  0.14422809;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKMNR2026", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }


                {
                    Options o{ { "n-resonances"_ok, "1"_ov } };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);

                    const complex<double> s_p = 4 * p["mass::pi^+"].evaluate() * p["mass::pi^+"].evaluate();

                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(s_p)),                                               -0.25748642,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(s_p)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.5)),                                               -0.04184462,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.5)),                                                0.58967422,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(1.0)),                                                0.36776120,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(1.0)),                                                0.92992027,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(complex<double>(0.5, 0.5))),                          0.24595382,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(complex<double>(0.5, 0.5))),                         -0.46087144,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(0.0)),                                                  2.19567550,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(0.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(0.5)),                                                  0.46958688,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(0.5)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0.5, 0.5))),                           -0.47140678,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0.5, 0.5))),                           -0.85748246,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(1.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(1.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(-1.0)),                                                 0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(-1.0)),                                                 0.00000000,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.0)),                                            -3.09366764,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.0)),                                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.5)),                                            -2.26285211,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.5)),                                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(complex<double>(0.5, 0.5))),                      -0.12448945,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(complex<double>(0.5, 0.5))),                       3.75965300,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(0,0.0)),                                    -3.09366764,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(0,0.0)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,0.0)),                                     2.19567550,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,0.0)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(0,0.5)),                                    -2.26285211,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(0,0.5)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,0.5)),                                    -0.66183918,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,0.5)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),              -2.41347800,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),               0.96009931,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(2,complex<double>(0.5, 0.5))),              -1.49375081,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(2,complex<double>(0.5, 0.5))),              -1.39113396,    eps);


                    const auto constrained_a = ff.constrained_a_fp_I1();
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[0],                                                 0.45544071,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[1],                                                -0.55225175,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[2],                                                -1.04508722,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[3],                                                -1.35538254,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.0)),                                                1.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.5)),                                                4.93318266,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.5)),                                                8.63581835,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(complex<double>(0.5, 0.5))),                         -0.80163832,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(complex<double>(0.5, 0.5))),                          1.58310240,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(1.0)),                                               -4.80021428,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(1.0)),                                               -0.75109616,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_21(complex<double>(0.0, 0.0))),                       0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_21(complex<double>(0.0, 0.0))),                       0.00000000,    eps);


                    //TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(complex<double>(0.0, 0.0))),                 -0.46788986,    eps);
                    //TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(complex<double>(0.0, 0.0))),                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(s_p)),                                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(s_p)),                                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(complex<double>(0.5, 0.5))),                  0.06802910,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(complex<double>(0.5, 0.5))),                  0.32270009,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(1.0)),                                        0.10012765,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(1.0)),                                        0.03951199,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_lenght_parameters()[0]),                       -7.79161435,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_lenght_parameters()[0]),                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_lenght_parameters()[1]),                      -15.04872684,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_lenght_parameters()[1]),                        0.00000000,    eps);

                    //needed for charged pion radius
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.0)),                                          -2.62154782,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.0)),                                           0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(s_p)),                                           0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(s_p)),                                           0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(complex<double>(0.5,0.5))),                     -4.06730590,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(complex<double>(0.5,0.5))),                     -7.25429346,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)),                                           5.86690477,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(1.0)),                                           8.65695683,    eps);
                }


                p["0->pipi::M_(+,1,1)@BHKMNR2026"]     =  1.465;
                p["0->pipi::Gamma_(+,1,1)@BHKMNR2026"] =  0.4;
                p["0->pipi::M_(+,1,2)@BHKMNR2026"]     =  1.72;
                p["0->pipi::Gamma_(+,1,2)@BHKMNR2026"] =  0.25;

                {
                    Options o{ { "n-resonances"_ok, "3"_ov } };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi22(complex<double>(2.10622500,-0.58600000))),            0.84434702,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi22(complex<double>(2.10622500,-0.58600000))),           -0.72846314,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi22(complex<double>(2.94277500,-0.43000000))),            0.86111963,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi22(complex<double>(2.94277500,-0.43000000))),           -0.59889900,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(0.0)),                                                  1.60479607,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(0.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0.5, 0.5))),                            3.98644434,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0.5, 0.5))),                           -1.57264136,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.0)),                                             2.43018174,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.0)),                                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(complex<double>(0.5, 0.5))),                       3.22772390,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(complex<double>(0.5, 0.5))),                     -12.46684166,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),              11.83372712,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),              -6.19220025,    eps);

                    const auto constrained_a2 = ff.constrained_a_fp_I1();
                    TEST_CHECK_NEARLY_EQUAL(constrained_a2[0],                                                 0.62313213,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a2[1],                                                -1.95879308,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a2[2],                                                 3.41500679,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a2[3],                                                -3.51275820,    eps);

                }
            }
        }
} parametric_BHKMNR2026_test;
