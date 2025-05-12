;;; Directory Local Variables            -*- no-byte-compile: t -*-
;;; For more information see (info "(emacs) Directory Variables")

((org-mode    . ((eval (lambda ()
                         (let ((project-root-dir (project-root (project-current))))
                           (setq-local python-shell-interpreter
                                       (concat project-root-dir "make.sh")
                                       python-shell-interpreter-args
                                       (concat "shell --share=/tmp -- python3 -i")
                                       python-shell-extra-pythonpaths
                                       (append
                                        (mapcar 'expand-file-name
                                                (list
                                                 (concat project-root-dir "src")
                                                 (concat project-root-dir "test")))
                                        python-shell-extra-pythonpaths)
                                       python-shell-prompt-detect-failure-warning
                                       nil))))))
 (python-mode . ((eval (lambda ()
                         (let ((project-root-dir (project-root (project-current))))
                           (setq-local python-shell-interpreter
                                       (concat project-root-dir "make.sh")
                                       python-shell-interpreter-args
                                       (concat "shell --share=/tmp -- python3 -i")
                                       python-shell-extra-pythonpaths
                                       (append
                                        (mapcar 'expand-file-name
                                                (list
                                                 (concat project-root-dir "src")
                                                 (concat project-root-dir "test")))
                                        python-shell-extra-pythonpaths)
                                       python-shell-prompt-detect-failure-warning
                                       nil)))))))
