GENCCODE_ASSEMBLY_TYPE=-a gcc
SO=so
SOBJ=so
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.58.2
COMPILE=/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/bin/x86_64-conda_cos6-linux-gnu-cc -D_FORTIFY_SOURCE=2 -O2 -D_REENTRANT  -DU_HAVE_ELF_H=1 -DU_HAVE_ATOMIC=1 -DU_HAVE_STRTOD_L=1  -DU_ATTRIBUTE_DEPRECATED= -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -pipe -I/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/include -std=c99 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -c
LIBFLAGS=-I/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/include -DPIC -fPIC
GENLIB=/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/bin/x86_64-conda_cos6-linux-gnu-cc -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -pipe -I/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/include -std=c99 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,-rpath,/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/lib -L/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/lib  -shared -Wl,-Bsymbolic
LDICUDTFLAGS=-nodefaultlibs -nostdlib
LD_SONAME=-Wl,-soname -Wl,
RPATH_FLAGS=
BIR_LDFLAGS=-Wl,-Bsymbolic
AR=/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/bin/x86_64-conda_cos6-linux-gnu-ar
ARFLAGS=r
RANLIB=/DATA/usr/n.klaassen/projects/SuRE_K562/analyses/SuRE_Pipeline/621cf0af/bin/x86_64-conda_cos6-linux-gnu-ranlib
INSTALL_CMD=/usr/bin/install -c
