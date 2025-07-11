(define-module (python-pyretechnics-package)
  #:use-module ((gnu packages base)            #:select (coreutils))
  #:use-module ((gnu packages certs)           #:select (nss-certs))
  #:use-module ((gnu packages check)           #:select (python-pytest
                                                         python-pytest-cov))
  #:use-module ((gnu packages emacs)           #:select (emacs-minimal))
  #:use-module ((gnu packages emacs-xyz)       #:select (emacs-htmlize))
  #:use-module ((gnu packages python)          #:select (python-wrapper))
  #:use-module ((gnu packages python-build)    #:select (python-setuptools-next
                                                         python-wheel
                                                         python-pypa-build
                                                         python-pyparsing))
  #:use-module ((gnu packages python-science)  #:select (python-distributed))
  #:use-module ((gnu packages python-web)      #:select (python-tornado-6))
  #:use-module ((gnu packages python-xyz)      #:select (python-numpy
                                                         python-rasterio
                                                         python-matplotlib
                                                         python-sortedcontainers
                                                         python-cython-3
                                                         python-twine
                                                         python-dask
                                                         python-fsspec))
  #:use-module ((gnu packages tls)             #:select (openssl))
  #:use-module ((gnu packages version-control) #:select (git))
  #:use-module ((guix build-system python)     #:select (pypi-uri))
  #:use-module ((guix build-system pyproject)  #:select (pyproject-build-system))
  #:use-module ((guix download)                #:select (url-fetch))
  #:use-module ((guix gexp)                    #:select (local-file))
  #:use-module ((guix git-download)            #:select (git-predicate))
  #:use-module ((guix licenses)                #:select (epl2.0 bsd-3 expat))
  #:use-module ((guix packages)                #:select (package
                                                         origin
                                                         base32
                                                         package-propagated-inputs
                                                         package-native-inputs))
  #:use-module ((guix utils)                   #:select (current-source-directory)))

;;============================================================
;; Utility Functions
;;============================================================

(define vcs-file?
  ;; Return true if the given file is under version control.
  (or (git-predicate (dirname (dirname (current-source-directory))))
      (const #t)))

;;============================================================
;; Dependency Packages
;;============================================================

(define-public python-snakeviz
  (package
   (name "python-snakeviz")
   (version "2.2.2")
   (source (origin
            (method url-fetch)
            (uri (pypi-uri "snakeviz" version))
            (sha256 (base32 "10vx7b0rn3gams2qk0dj5xkbiavh8x13hykm2kzk581lirpqq0h8"))))
   (build-system pyproject-build-system)
   (arguments `(#:tests? #f))
   (propagated-inputs (list python-tornado-6))
   (native-inputs (list python-setuptools-next python-wheel))
   (home-page "https://jiffyclub.github.io/snakeviz/")
   (synopsis "A web-based viewer for Python profiler output")
   (description "This package provides a web-based viewer for Python profiler output.")
   (license bsd-3)))

(define-public python-pkginfo
  (package
   (name "python-pkginfo")
   (version "1.12.1.2")
   (source (origin
            (method url-fetch)
            (uri (pypi-uri "pkginfo" version))
            (sha256 (base32 "0ysbl0rhy8chvx58zi5qb61jli76dcyblklnc0118vy39a15gnaw"))))
   (build-system pyproject-build-system)
   (native-inputs (list python-pytest python-pytest-cov python-setuptools-next python-wheel))
   (home-page "https://code.launchpad.net/~tseaver/pkginfo/trunk")
   (synopsis "Query metadata from sdists / bdists / installed packages.")
   (description "Query metadata from sdists / bdists / installed packages.")
   (license expat)))

(define-public python-rasterio-next
  (package
   (inherit python-rasterio)
   (version "1.4.3")
   (source (origin
            (method url-fetch)
            (uri (pypi-uri "rasterio" version))
            (sha256 (base32 "02immzrf8g3ms75nrld4n1jp42www3s1r2n7nan9swy4qzdha7r0"))))
   (propagated-inputs (cons python-pyparsing
                            (map cadr
                                 (assoc-remove!
                                  (assoc-remove!
                                   (package-propagated-inputs python-rasterio)
                                   "python-snuggs")
                                  "python-setuptools"))))
   (native-inputs (cons* python-fsspec
                         python-numpy
                         python-setuptools-next
                         python-cython-3
                         (map cadr
                              (assoc-remove!
                               (package-native-inputs python-rasterio)
                               "python-cython"))))))

;;============================================================
;; The Pyretechnics Package
;;============================================================

(define-public python-pyretechnics
  (package
   (name "python-pyretechnics")
   (version "2025.6.13")
   (source (local-file "../.."
                       "pyretechnics-checkout"
                       #:recursive? #t
                       #:select?    vcs-file?))
   (build-system pyproject-build-system) ; includes python-sans-pip-wrapper
   ;; See /run/current-system/profile/share/guile/site/3.0/guix/build/pyproject-build-system.scm for more info.
   (arguments `(#:configure-flags '() ; passed to setuptools.build_meta.build_wheel(wheel_dir, config_settings=$(configure-flags))
                #:tests?          #f
                #:test-backend    'pytest
                #:test-flags      '())) ; passed to pytest -vv $(test-flags)
   (native-inputs (list
                   ;; Shell utilities
                   coreutils
                   ;; SSL/TLS
                   openssl
                   nss-certs
                   ;; Literate programming
                   emacs-minimal
                   emacs-htmlize
                   git ; used by emacs
                   ;; Python3
                   python-wrapper
                   ;; Build/Test/Profile/Deploy utilities
                   python-setuptools-next
                   python-wheel
                   python-pypa-build
                   python-pytest
                   python-snakeviz
                   python-twine
                   python-pkginfo
                   ;; Notebook programming examples
                   python-matplotlib
                   python-dask
                   python-distributed))
   (propagated-inputs (list
                       ;; Runtime dependencies
                       python-numpy
                       python-rasterio-next
                       python-sortedcontainers
                       python-cython-3))
   (synopsis "A Python library for simulating fire behavior in a variety of ways.")
   (description "A Python library for simulating fire behavior in a variety of ways.")
   (home-page "https://github.com/pyregence/pyretechnics/")
   (license epl2.0)))

;; Return the package for guix commands that expect it
python-pyretechnics
