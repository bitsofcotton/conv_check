CXX=	clang++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -mtune=native
#CXXFLAGS+=	-g2 -O2 -Wall
#CXXFLAGS+=	-fopenmp -lgomp
#CXXFLAGS+=	-pg
LDFLAGS=	-lc++
#LDFLAGS+=	-L/usr/local/lib -lmpfr -lgmp
#LDFLAGS+=	-L/usr/local/lib -lqd -lc++

# extra easy problems.
#CXXFLAGS+=	-DACC_DOUBLE
CXXFLAGS+=	-DACC_LDOUBLE

# for normal use.
#CXXFLAGS+=	-DACC_QD_DDOUBLE
#CXXFLAGS+=	-DACC_QD_QDOUBLE
#CXXFLAGS+=	-DACC_GMP=256
#CXXFLAGS+=	-DACC_GMP=512

# with integer only (too slow)
#CXXFLAGS+=	-DACC_NO_FLOAT

# Without eigen. do not use this because it costs dramatically long time.
# N.B. not using with cpu implemented float, it costs slight long time.
CXXFLAGS+=	-D_WITHOUT_EIGEN_

CLEANFILES+=	konbu

clean:
	@rm -rf ${CLEANFILES}

