;;; Directory Local Variables            -*- no-byte-compile: t -*-
;;; For more information see (info "(emacs) Directory Variables")

((org-mode    . ((eval (lambda ()
                         (let ((project-root-dir (project-root (project-current))))
                           (setq-local python-shell-extra-pythonpaths
                                       (append
                                        (mapcar 'expand-file-name
                                                (list
                                                 (concat project-root-dir "src")
                                                 (concat project-root-dir "test")))
                                        python-shell-extra-pythonpaths)))))))
 (python-mode . ((eval (lambda ()
                         (let ((project-root-dir (project-root (project-current))))
                           (setq-local python-shell-extra-pythonpaths
                                       (append
                                        (mapcar 'expand-file-name
                                                (list
                                                 (concat project-root-dir "src")
                                                 (concat project-root-dir "test")))
                                        python-shell-extra-pythonpaths))))))))
