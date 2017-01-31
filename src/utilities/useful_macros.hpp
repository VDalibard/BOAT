#ifndef USEFUL_MACROS_HPP_INCLUDED
#define USEFUL_MACROS_HPP_INCLUDED

#include "boost/preprocessor/arithmetic/sub.hpp"
#include <iostream>

#define PP_NARG_M1(...) BOOST_PP_SUB(PP_NARG(__VA_ARGS__), 1)

#define PP_NARG(...) PP_NARG_(_0, ##__VA_ARGS__, PP_RSEQ_N())
#define PP_NARG_(...) PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, \
                 _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26,  \
                 _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38,  \
                 _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50,  \
                 _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62,  \
                 _63, _64, N, ...)                                            \
  N
#define PP_RSEQ_N()                                                           \
  63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, \
      44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, \
      26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9,  \
      8, 7, 6, 5, 4, 3, 2, 1, 0

#define _PP_0(_1, ...) _1             // (a,b,c,d) => a
#define _PP_X(_1, ...) (__VA_ARGS__)  // (a,b,c,d) => (b,c,d)

#define PP_MAP_0(f, a)
#define PP_MAP_1(f, a) f(_PP_0 a) PP_MAP_0(f, _PP_X a)
#define PP_MAP_2(f, a) f(_PP_0 a) PP_MAP_1(f, _PP_X a)
#define PP_MAP_3(f, a) f(_PP_0 a) PP_MAP_2(f, _PP_X a)
#define PP_MAP_4(f, a) f(_PP_0 a) PP_MAP_3(f, _PP_X a)
#define PP_MAP_5(f, a) f(_PP_0 a) PP_MAP_4(f, _PP_X a)
#define PP_MAP_6(f, a) f(_PP_0 a) PP_MAP_5(f, _PP_X a)
#define PP_MAP_7(f, a) f(_PP_0 a) PP_MAP_6(f, _PP_X a)
#define PP_MAP_8(f, a) f(_PP_0 a) PP_MAP_7(f, _PP_X a)
#define PP_MAP_9(f, a) f(_PP_0 a) PP_MAP_8(f, _PP_X a)
#define PP_MAP_10(f, a) f(_PP_0 a) PP_MAP_9(f, _PP_X a)
#define PP_MAP_11(f, a) f(_PP_0 a) PP_MAP_10(f, _PP_X a)
#define PP_MAP_12(f, a) f(_PP_0 a) PP_MAP_11(f, _PP_X a)
#define PP_MAP_13(f, a) f(_PP_0 a) PP_MAP_12(f, _PP_X a)
#define PP_MAP_14(f, a) f(_PP_0 a) PP_MAP_13(f, _PP_X a)
#define PP_MAP_15(f, a) f(_PP_0 a) PP_MAP_14(f, _PP_X a)
#define PP_MAP_16(f, a) f(_PP_0 a) PP_MAP_15(f, _PP_X a)
#define PP_MAP_17(f, a) f(_PP_0 a) PP_MAP_16(f, _PP_X a)
#define PP_MAP_18(f, a) f(_PP_0 a) PP_MAP_17(f, _PP_X a)
#define PP_MAP_19(f, a) f(_PP_0 a) PP_MAP_18(f, _PP_X a)
#define PP_MAP_20(f, a) f(_PP_0 a) PP_MAP_19(f, _PP_X a)
#define PP_MAP_21(f, a) f(_PP_0 a) PP_MAP_20(f, _PP_X a)
#define PP_MAP_22(f, a) f(_PP_0 a) PP_MAP_21(f, _PP_X a)
#define PP_MAP_23(f, a) f(_PP_0 a) PP_MAP_22(f, _PP_X a)
#define PP_MAP_24(f, a) f(_PP_0 a) PP_MAP_23(f, _PP_X a)
#define PP_MAP_25(f, a) f(_PP_0 a) PP_MAP_24(f, _PP_X a)
#define PP_MAP_26(f, a) f(_PP_0 a) PP_MAP_25(f, _PP_X a)
#define PP_MAP_27(f, a) f(_PP_0 a) PP_MAP_26(f, _PP_X a)
#define PP_MAP_28(f, a) f(_PP_0 a) PP_MAP_27(f, _PP_X a)
#define PP_MAP_29(f, a) f(_PP_0 a) PP_MAP_28(f, _PP_X a)
#define PP_MAP_30(f, a) f(_PP_0 a) PP_MAP_29(f, _PP_X a)
#define PP_MAP_31(f, a) f(_PP_0 a) PP_MAP_30(f, _PP_X a)
#define PP_MAP_32(f, a) f(_PP_0 a) PP_MAP_31(f, _PP_X a)
#define PP_MAP_33(f, a) f(_PP_0 a) PP_MAP_32(f, _PP_X a)
#define PP_MAP_34(f, a) f(_PP_0 a) PP_MAP_33(f, _PP_X a)
#define PP_MAP_35(f, a) f(_PP_0 a) PP_MAP_34(f, _PP_X a)
#define PP_MAP_36(f, a) f(_PP_0 a) PP_MAP_35(f, _PP_X a)
#define PP_MAP_37(f, a) f(_PP_0 a) PP_MAP_36(f, _PP_X a)
#define PP_MAP_38(f, a) f(_PP_0 a) PP_MAP_37(f, _PP_X a)
#define PP_MAP_39(f, a) f(_PP_0 a) PP_MAP_38(f, _PP_X a)
#define PP_MAP_40(f, a) f(_PP_0 a) PP_MAP_39(f, _PP_X a)
#define PP_MAP_41(f, a) f(_PP_0 a) PP_MAP_40(f, _PP_X a)
#define PP_MAP_42(f, a) f(_PP_0 a) PP_MAP_41(f, _PP_X a)
#define PP_MAP_43(f, a) f(_PP_0 a) PP_MAP_42(f, _PP_X a)
#define PP_MAP_44(f, a) f(_PP_0 a) PP_MAP_43(f, _PP_X a)
#define PP_MAP_45(f, a) f(_PP_0 a) PP_MAP_44(f, _PP_X a)
#define PP_MAP_46(f, a) f(_PP_0 a) PP_MAP_45(f, _PP_X a)
#define PP_MAP_47(f, a) f(_PP_0 a) PP_MAP_46(f, _PP_X a)
#define PP_MAP_48(f, a) f(_PP_0 a) PP_MAP_47(f, _PP_X a)
#define PP_MAP_49(f, a) f(_PP_0 a) PP_MAP_48(f, _PP_X a)
#define PP_MAP_50(f, a) f(_PP_0 a) PP_MAP_49(f, _PP_X a)
#define PP_MAP_51(f, a) f(_PP_0 a) PP_MAP_50(f, _PP_X a)
#define PP_MAP_52(f, a) f(_PP_0 a) PP_MAP_51(f, _PP_X a)
#define PP_MAP_53(f, a) f(_PP_0 a) PP_MAP_52(f, _PP_X a)
#define PP_MAP_54(f, a) f(_PP_0 a) PP_MAP_53(f, _PP_X a)
#define PP_MAP_55(f, a) f(_PP_0 a) PP_MAP_54(f, _PP_X a)
#define PP_MAP_56(f, a) f(_PP_0 a) PP_MAP_55(f, _PP_X a)
#define PP_MAP_57(f, a) f(_PP_0 a) PP_MAP_56(f, _PP_X a)
#define PP_MAP_58(f, a) f(_PP_0 a) PP_MAP_57(f, _PP_X a)
#define PP_MAP_59(f, a) f(_PP_0 a) PP_MAP_58(f, _PP_X a)
#define PP_MAP_60(f, a) f(_PP_0 a) PP_MAP_59(f, _PP_X a)
#define PP_MAP_61(f, a) f(_PP_0 a) PP_MAP_60(f, _PP_X a)
#define PP_MAP_62(f, a) f(_PP_0 a) PP_MAP_61(f, _PP_X a)
#define PP_MAP_63(f, a) f(_PP_0 a) PP_MAP_62(f, _PP_X a)
#define PP_MAP_64(f, a) f(_PP_0 a) PP_MAP_63(f, _PP_X a)

#define PP_MAP(f, ...) \
  CAT(PP_MAP_, PP_NARG_M1(_0, ##__VA_ARGS__))(f, (__VA_ARGS__))

#define CAT(a, ...) PRIMITIVE_CAT(a, __VA_ARGS__)
#define PRIMITIVE_CAT(a, ...) a##__VA_ARGS__

#define PR(...) std::cout PP_MAP(PR_, ##__VA_ARGS__) << "\n"
#define PR_(A) << QUOTE(A) << " = " << A << ",   "
#define QUOTE(A) #A
#define PRS(S) std::cout << #S << "\n"
#define SIZE(A) int((A).size())
#define BLOCK()                 \
  {                             \
    int ______a_______;         \
    std::cin >> ______a_______; \
  }

#define RST "\x1B[0m"
#define KRED "\x1B[31m"
#define KGRN "\x1B[32m"
#define KYEL "\x1B[33m"
#define KBLU "\x1B[34m"
#define KMAG "\x1B[35m"
#define KCYN "\x1B[36m"
#define KWHT "\x1B[37m"
#define BOLD "\x1B[1m"

#define AVG_PROP(ENGINE, PROPERTY) ENGINE.average([&](const auto& m){return m.PROPERTY;})


#endif  // USEFUL_MACROS_HPP_INCLUDED
