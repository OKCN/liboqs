AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I config

# DOXYGEN SUPPORT
include aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg

CLANGFORMAT ?= clang-format

SUBDIRS = ${SRCDIR}

BUILT_SOURCES = links
lib_LTLIBRARIES = liboqs.la
liboqs_la_SOURCES =
liboqs_la_LIBADD = src/common/libcommon.la src/kex/libkex.la src/crypto/rand/librand.la src/crypto/aes/libaes.la
liboqs_la_LIBADD += src/crypto/rand_urandom_aesctr/librandaesctr.la src/crypto/sha3/libsha3.la
liboqs_la_LIBADD += src/crypto/rand_urandom_chacha20/librandchacha20.la
liboqs_la_LIBADD += src/kex_rlwe_bcns15/libbcns15.la src/kex_rlwe_newhope/libnewhope.la 
liboqs_la_LIBADD += src/kex_rlwe_newhope_simple/libnewhopesimple.la src/kex_rlwe_zarzar/libzarzar.la src/kex_rlwe_E8/libE8.la
liboqs_la_LIBADD += src/kex_lwe_frodo/libfrodo.la src/kex_rlwe_msrln16/libmsrln16.la
liboqs_la_LIBADD += src/kex_lwe_okcn/libokcn.la src/kex_lwr_okcn/liblwrokcn.la
liboqs_la_LIBADD += src/kex_sidh_cln16/libcln16.la

if USE_MCBITS
liboqs_la_LIBADD += src/kex_code_mcbits/libmcbits.la
endif

if USE_NTRU
liboqs_la_LIBADD += src/kex_ntru/libntru.la
liboqs_la_LIBADD += external/NTRUEncrypt-master/.libs/libntruencrypt.la
endif

if USE_SIDH_IQC
liboqs_la_LIBADD += src/kex_sidh_iqc_ref/libsidhiqc.la
endif

noinst_bin_PROGRAMS = test_rand test_kex test_aes
noinst_bindir=$(prefix)/tests
test_kex_LDADD =  liboqs.la -lm
test_kex_SOURCES = src/kex/test_kex.c
test_kex_CPPFLAGS = -I./include
test_kex_CPPFLAGS += $(AM_CPPFLAGS)
if USE_OPENSSL
test_kex_LDADD += -lcrypto
endif
if USE_MCBITS
test_kex_LDADD += -lsodium
endif
if USE_SIDH_IQC
test_kex_LDADD += -lgmp
endif

test_aes_LDADD =  liboqs.la -lm
test_aes_SOURCES = src/crypto/aes/test_aes.c
test_aes_CPPFLAGS = -I./include
test_aes_CPPFLAGS += $(AM_CPPFLAGS)


if USE_OPENSSL
test_aes_LDADD += -lcrypto
endif

test_rand_SOURCES = src/crypto/rand/test_rand.c

test_rand_CPPFLAGS = -Iinclude -Isrc/crypto/rand_urandom_aesctr/
test_rand_CPPFLAGS += $(AM_CPPFLAGS)
test_rand_LDADD = liboqs.la
if USE_OPENSSL
test_rand_LDADD += -lcrypto
else
if USE_AES_NI
test_rand_CPPFLAGS += -maes -msse2
endif
endif


test: clean-tests
	make
if USE_NTRU
if ON_DARWIN
	install_name_tool -change `otool -L ./test_kex | grep "ntru" | awk -F ' ' '{ print $$1 }'` external/NTRUEncrypt-master/.libs/libntruencrypt.0.dylib ./test_kex
	install_name_tool -change `otool -L ./test_aes | grep "ntru" | awk -F ' ' '{ print $$1 }'`  external/NTRUEncrypt-master/.libs/libntruencrypt.0.dylib ./test_aes
	install_name_tool -change `otool -L ./test_rand | grep "ntru" | awk -F ' ' '{ print $$1 }'` external/NTRUEncrypt-master/.libs/libntruencrypt.0.dylib ./test_rand
endif
endif
	./test_kex
	./test_rand
	./test_aes

links:
	$(MKDIR_P) include/oqs
	$(LN_S) -f ../../src/common/common.h include/oqs
	$(LN_S) -f ../../src/crypto/aes/aes.h include/oqs
	$(LN_S) -f ../../src/crypto/sha3/sha3.h include/oqs
	$(LN_S) -f ../../src/kex/kex.h include/oqs
	$(LN_S) -f ../../src/crypto/rand/rand.h include/oqs
	$(LN_S) -f ../../src/crypto/rand_urandom_chacha20/rand_urandom_chacha20.h include/oqs
	$(LN_S) -f ../../src/crypto/rand_urandom_aesctr/rand_urandom_aesctr.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_newhope/kex_rlwe_newhope.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_newhope_simple/kex_rlwe_newhope_simple.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_zarzar/kex_rlwe_zarzar.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_E8/kex_rlwe_E8.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_msrln16/kex_rlwe_msrln16.h include/oqs
	$(LN_S) -f ../../src/kex_lwe_frodo/kex_lwe_frodo.h include/oqs
	$(LN_S) -f ../../src/kex_lwe_okcn/kex_lwe_okcn.h include/oqs
	$(LN_S) -f ../../src/kex_lwr_okcn/kex_lwr_okcn.h include/oqs
	$(LN_S) -f ../../src/kex_rlwe_bcns15/kex_rlwe_bcns15.h include/oqs
	$(LN_S) -f ../../src/kex_sidh_cln16/kex_sidh_cln16.h include/oqs
	$(LN_S) -f .libs/liboqs.a
if USE_MCBITS
	$(LN_S) -f ../../src/kex_code_mcbits/kex_code_mcbits.h include/oqs
endif
if USE_NTRU
	$(LN_S) -f ../../src/kex_ntru/kex_ntru.h include/oqs
	test -d "external/NTRUEncrypt-master" || ./download-and-build-ntru.sh
endif
if USE_SIDH_IQC
	$(LN_S) -f ../../src/kex_sidh_iqc_ref/kex_sidh_iqc_ref.h include/oqs
	$(LN_S) -f src/kex_sidh_iqc_ref/sample_params
endif

clean-local:
	rm -f liboqs.a

clean-tests:
	rm -f test_kex test_rand test_aes

prettyprint:
	find src -name '*.c' -o -name '*.h' | xargs $(CLANGFORMAT) -style=file -i

docs: links
	doxygen
