#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Load pyretechnics.org and tangle its code blocks to disk
;;==========================================================

(require 'org)
(require 'python-mode)

(find-file "pyretechnics.org")
(org-babel-tangle)

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
