/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CONFIG_H
#define GMATELASTOPLASTICFINITESTRAINSIMO_CONFIG_H

#ifdef GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT

    #define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(expr) \
        GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT_IMPL(expr, __FILE__, __LINE__)

    #define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }

#else

    #define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(expr)

#endif

#endif
