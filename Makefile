# gcc.
CXX=	g++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-O3 -mtune=native
#CXXFLAGS+=	-g2 -O0 -Wall
#CXXFLAGS+=	-fopenmp -pthread
#CXXFLAGS+=	-pg
LDFLAGS=	-lstdc++
#LDFLAGS+=	-L/usr/local/lib -lmpfr -lgmp

# test for tiny dimensions.
#CXXFLAGS+=	-DACC_FLOAT

# normal use.
CXXFLAGS+=	-DACC_DOUBLE
#CXXFLAGS+=	-DACC_LDOUBLE

# for hard problems including netlib's small-sized problems with 0<=x condition.
#CXXFLAGS+=	-DACC_GMP=256
#CXXFLAGS+=	-DACC_GMP=512

# Without eigen, do not use this because it costs dramatically long time.
CXXFLAGS+=	-DWITHOUT_EIGEN

CLEANFILES+=	konbu

clean:
	@rm -rf ${CLEANFILES}

