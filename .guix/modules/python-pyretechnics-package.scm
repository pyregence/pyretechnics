(define-module (python-pyretechnics-package)
  #:use-module ((gnu packages base)            #:select (coreutils which))
  #:use-module ((gnu packages bash)            #:select (bash))
  #:use-module ((gnu packages check)           #:select (python-pytest))
  #:use-module ((gnu packages emacs)           #:select (emacs-minimal))
  #:use-module ((gnu packages emacs-xyz)       #:select (emacs-htmlize))
  #:use-module ((gnu packages geo)             #:select (gdal))
  #:use-module ((gnu packages less)            #:select (less))
  #:use-module ((gnu packages python)          #:select (python-wrapper))
  #:use-module ((gnu packages python-build)    #:select (python-hatchling))
  #:use-module ((gnu packages python-xyz)      #:select (python-numpy python-rasterio python-matplotlib python-twine))
  #:use-module ((gnu packages ssh)             #:select (openssh))
  #:use-module ((gnu packages version-control) #:select (git))
  #:use-module ((guix build-system pyproject)  #:select (pyproject-build-system))
  #:use-module ((guix gexp)                    #:select (local-file))
  #:use-module ((guix git-download)            #:select (git-predicate))
  #:use-module ((guix licenses)                #:select (epl2.0))
  #:use-module ((guix packages)                #:select (package))
  #:use-module ((guix utils)                   #:select (current-source-directory))
  #:use-module ((lambdatronic python-packages) #:select (python-snakeviz python-cython-3)))

(define vcs-file?
  ;; Return true if the given file is under version control.
  (or (git-predicate (dirname (dirname (current-source-directory))))
      (const #t)))

(define-public python-pyretechnics
  (package
   (name "python-pyretechnics")
   (version "2024.11.7")
   (source (local-file "../.."
                       "pyretechnics-checkout"
                       #:recursive? #t
                       #:select?    vcs-file?))
   (build-system pyproject-build-system) ; includes python-toolchain
   ;; See file:/run/current-system/profile/share/guile/site/3.0/guix/build/pyproject-build-system.scm for more info.
   (arguments `(#:configure-flags '()
                #:test-backend    'pytest
                #:test-flags      '()))
   (native-inputs (list
                   ;; Shell utilities
                   bash
                   coreutils
                   which
                   less
                   ;; Version control
                   git
                   openssh
                   ;; GIS utilities
                   gdal
                   ;; Build tools
                   emacs-minimal
                   emacs-htmlize
                   python-hatchling
                   python-pytest
                   python-twine
                   python-snakeviz))
   (propagated-inputs (list
                       ;; Python3
                       python-wrapper
                       ;; Python dependency libraries
                       python-numpy
                       python-rasterio
                       python-matplotlib
                       python-cython-3))
   (synopsis "A Python library for simulating fire behavior in a variety of ways.")
   (description "A Python library for simulating fire behavior in a variety of ways.")
   (home-page "https://github.com/pyregence/pyretechnics/")
   (license epl2.0)))

;; Return the package for guix commands that expect it
python-pyretechnics
