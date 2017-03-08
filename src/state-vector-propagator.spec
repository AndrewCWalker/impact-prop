Name:           state-vector-propagator
Version:        0.1plus
Release:        1%{?dist}
Summary:        state vector propagator library from the impact project

License:        GPL (FIXME: check this)
URL:            http://impact.lanl.gov/
Source0:        %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

#BuildRequires:  tex

%description
The state vector propagator library from the impact project.  Currently 
built without using cspice.

%prep
%setup -q

%build
cd src
make %{?_smp_mflags} -f MakefileDIORAMA

%install
cd src
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/%{_bindir}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}
mkdir -p $RPM_BUILD_ROOT/%{_includedir}
cp propagator $RPM_BUILD_ROOT/%{_bindir}
cp libpropagator.so $RPM_BUILD_ROOT/%{_libdir}
cp defs.h structs.h test.h acceleration.h assoclegendre.h atmprop.h Cd.h emu.h interpolation.h io.h misc.h msisinputs.h nrlmsise-00.h propagator.h rk.h transformation.h $RPM_BUILD_ROOT/%{_includedir}


%files
%doc README
%{_includedir}/*.h
%{_bindir}/*
%{_libdir}/*

%changelog
* Wed Jan 21 2015 Mark Galassi
- first packaging
