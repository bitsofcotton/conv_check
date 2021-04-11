CXX=	clang++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -gfull -mtune=native
#CXXFLAGS+=	-fopenmp -L/usr/local/lib -lgomp
#CXXFLAGS+=	-pg
LDFLAGS=	-lc++

#CXXFLAGS+=     -D_FLOAT_BITS_=32
#CXXFLAGS+=     -D_FLOAT_BITS_=256

CLEANFILES+=	konbu

clean:
	@rm -rf ${CLEANFILES}

