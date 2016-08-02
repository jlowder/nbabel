(ql:quickload :cl-ppcre :silent t)
(use-package :cl-ppcre)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; some general arbitrary-length vector utility functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vlen^2 (r &optional (acc 0))
  "length of a vector, squared"
  (cond ((null r) acc)
        ((numberp r) (+ acc (* r r)))
        (t (vlen^2 (cdr r) (+ acc (* (car r) (car r)))))))

(defun vlen (r)
  "length of a vector"
  (sqrt (vlen^2 r)))

(defun s*v (s v)
  "scalar times vector"
  (loop for i in v
     collect (* s i)))

(defun v-v (&rest r)
  "vector subtraction, first vector minus all remaining vectors passed in"
  (flet ((rec (p1 p2)
           (loop for a in p1
              for b in p2
              collect (- a b))))
    (reduce #'rec r)))

(defun v+v (&rest r)
  "vector addition, sum of all vectors passed in"
  (flet ((rec (v1 v2)
           (loop for x in v1
              for y in v2
              collect (+ x y))))
    (reduce #'rec r)))

(defun dot (a b &optional (acc 0))
  "dot product of two vectors"
  (if (or (null a) (null b))
      acc
      (dot (cdr a) (cdr b) (+ acc (* (car a) (car b))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; beginning of nbabel-specific code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun force (ri rj mj)
  "calculate the force between two particles"
  (let* ((p (v-v rj ri))
         (c (* mj (expt (vlen^2 p) -1.5))))
    (s*v c p)))

(defun file->strlist (filename)
  "read a file into a list of strings"
  (with-open-file (stream filename)
    (loop for line = (read-line stream nil)
       while line
       collect line)))

(defun ingest (fname)
  "parse an input file. return a list of positions, velocities, and masses (each as lists)."
  (let ((p1 (remove-if (lambda (x) (< (length x) 3)) (mapcar (lambda (x) (split " +" x)) (file->strlist fname)))))
    (loop for line in p1
       as parsed = (mapcar (lambda (x) (with-input-from-string (l x) (read l))) line)
       collecting (nth 1 parsed) into masses
       collecting (subseq parsed 2 5) into x-data
       collecting (subseq parsed 5) into v-data
       finally (return (list x-data v-data masses)))))

(defun lineout (time ke pe total err)
  (format t "~A~24T~A~48T~A~72T~A~96T~A~%" time ke pe total err))

(defun summary (time ke pe e0)
  (let ((te (- ke pe)))
    (lineout time ke pe te (/ (- te e0) e0))))

(defun main ()
  (destructuring-bind (xdata vdata masses)
      (ingest (or (car (uiop:command-line-arguments)) "input16"))
    (flet ((ke ()
             "kinetic energy calculation"
             (loop for v in vdata
                for m in masses
                sum (* 1/2 m (dot v v))))
           (pe ()
             "potential energy calculation"
             (let ((n (length xdata)))
               (loop for i from 0 to (1- n) sum
                    (loop for j from (1+ i) to (1- n)
                       sum (/ (* (nth i masses) (nth j masses)) (vlen (v-v (nth j xdata) (nth i xdata))))))))
           (leapfrog (f)
             "create a closure representing the leapfrog algorithm. `F` is a function that calculates a particle's instantaneous acceleration."
             (lambda (i dt)
               (let ((xi (nth i xdata))
                     (vi (nth i vdata)))
                 (let* ((ai (funcall f i xi))
                        (xi+1 (v+v xi (s*v dt vi) (s*v (* 1/2 dt dt) ai)))
                        (ai+1 (funcall f i xi+1)))
                   (values
                    xi+1
                    (v+v vi (s*v (* 1/2 dt) (v+v ai ai+1))))))))
           (Fi (i p)
             "calculate the total force acting on particle i, at position p"
             (loop for xi in xdata
                as j = 0 then (incf j)
                as mj = (nth j masses)
                as fi = (if (not (= i j))
                            (force p xi mj)
                            (list 0 0 0))
                as acc = fi then (v+v acc fi)
                finally (return acc))))
      (let ((integrator (leapfrog #'Fi))
            (dt 1d-3)
            (brk "=======================")
            (e0 (- (ke) (pe)))
            (tend (with-input-from-string (l (or (second (uiop:command-line-arguments)) "1")) (read l))))
        (lineout "Time" "Kinetic Energy" "Potential Energy" "Total Energy" "Energy Error")
        (lineout brk brk brk brk brk)
        (summary 0d0 (ke) (pe) e0)
        (loop for time from 0 to tend by dt
           as inc = 0 then (1+ inc)
           do (loop for i from 0 to (1- (length xdata))
                 do (multiple-value-bind (nx nv) (funcall integrator i dt)
                      (setf (nth i xdata) nx)
                      (setf (nth i vdata) nv)))
             (when (zerop (mod inc 10))
               (summary time (ke) (pe) e0)))))))
  
