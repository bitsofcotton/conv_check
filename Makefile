CXX=	clang++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -gfull -mtune=native
#CXXFLAGS+=	-Ofast -gfull -mno-sse -mno-sse2 -mno-sse3 -mno-3dnow -mno-mmx
CXXFLAGS+=	-fopenmp -L/usr/local/lib -lomp
#CXXFLAGS+=	-pg
LDFLAGS=	-lc++

#CXXFLAGS+=     -D_FLOAT_BITS_=32
#CXXFLAGS+=     -D_FLOAT_BITS_=64
#CXXFLAGS+=     -D_FLOAT_BITS_=128
#CXXFLAGS+=     -D_FLOAT_BITS_=256

CLEANFILES+=	konbu

clean:
	@rm -rf ${CLEANFILES}

