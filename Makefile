# gcc.
CXX=	g++
LD=	${CXX}

# compiler flags.
CXXFLAGS=	-I/usr/local/include
CXXFLAGS+=	-O3 -mtune=core2
#CXXFLAGS+=	-g2 -O0 -Wall
#CXXFLAGS+=	-fno-inline
#CXXFLAGS+=	-fopenmp -pthread
LDFLAGS=	-lstdc++
#LDFLAGS+=	-L/usr/local/lib -lmpfr -lgmp

# test for tiny dimensions.
#CXXFLAGS+=	-DACC_GMP -DBITS=14
#CXXFLAGS+=	-DACC_GMP -DBITS=15
#CXXFLAGS+=	-DACC_GMP -DBITS=16
#CXXFLAGS+=	-DACC_GMP -DBITS=17
#CXXFLAGS+=	-DACC_FLOAT

# normal use.
#CXXFLAGS+=	-DACC_GMP -DBITS=24
#CXXFLAGS+=	-DACC_GMP -DBITS=32
CXXFLAGS+=	-DACC_DDOUBLE
#CXXFLAGS+=	-DACC_LLDOUBLE

# for hard problems like bandm or degen and so on.
#CXXFLAGS+=	-DACC_GMP -DBITS=64
#CXXFLAGS+=	-DACC_GMP -DBITS=128
#CXXFLAGS+=	-DACC_GMP -DBITS=160
#CXXFLAGS+=	-DACC_GMP -DBITS=256
#CXXFLAGS+=	-DACC_GMP -DBITS=512

clean:
	@rm -rf ${CLEANFILES}

