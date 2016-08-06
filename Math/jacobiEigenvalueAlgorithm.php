<?php namespace PCA\Math;

/**
 * PCA\Math\jacobiEigenvalueAlgorithm
 *
 * Applications for real symmetric matrices only.
 * In other cases the above implementation will never terminate.
 *
 * @author laxip
 *
 * $example_matrix = array(array(  4, -30,   60,  -35),
 *                         array(-30, 300, -675,  420),
 *                         array( 60,-675, 1620,-1050),
 *                         array(-35, 420,-1050,  700));
 *
 */
class jacobiEigenvalueAlgorithm
{
    /**
     * @var array matrix
     */
    protected $S;

    /**
     * @var int|null size of matrix
     */
    protected $n;

    /**
     * matrix which contains coords for eigenvectors
     * @var array
     */
    protected $_E;

    /**
     * the corresponding eigenvectors
     * @var array
     */
    protected $V;

    /**
     * vector e which contains the eigenvalues
     * @var
     */
    protected $e;

    /**
     * @var array
     */
    protected $changed = array();

    /**
     * @var
     */
    protected $state;

    /**
     * @var
     */
    protected $iter;

    /**
     * jacobiEigenvalueAlgorithm constructor.
     * @param array $_S
     * @param null $_n
     */
    function __construct(array $_S, $_n = null)
    {
        $this->S = $_S;
        $this->n = (is_null($_n)) ? count($this->S) : $_n;

        /*
         * Add row on top and for each row add 0 at the beginning to
         * produce matrix which can be indexed by 1.
         */
        for ($i = 0; $i < $this->n; $i++) {
            array_unshift($this->S[$i], 0);
        }
        array_unshift($this->S, array(0));

        $this->calculate();
    }

    /**
     * start Jacobi algorithm
     */
    protected final function calculate()
    {
        /*
         * init e, E, and arrays ind, changed
         */
        $this->I();
        $this->state = $this->n;

        /*
         * fix for association
         */
        $ind[0] = $this->e[0] = $this->changed[0] = 0;

        for ($k = 1; $k <= $this->n; $k++) {
            $ind[$k] = $this->maxind($k);
            $this->e[$k] = $this->S[$k][$k];
            $this->changed[$k] = true;
        }


        $this->iter = 0;


        /*
         * next rotation
         */
        while ($this->state != 0) {

            /*
             * find index (k,l) of pivot p
             */
            $m = 1;

            for ($k = 2; $k <= $this->n - 1; $k++) {
                if (abs($this->S[$k][$ind[$k]]) >
                    abs($this->S[$m][$ind[$m]])
                ) {
                    $m = $k;
                }
            }

            $k = $m;
            $l = $ind[$m];
            $p = $this->S[$k][$l];
            $pp = $p * $p;

            /*
             * calculate c = cos φ, s = sin φ
             */
            $y = ($this->e[$l] - $this->e[$k]) / 2;
            $d = abs($y) + sqrt($pp + $y * $y);
            $r = sqrt($pp + $d * $d);

            $this->cos = $d / $r;
            $this->sin = $p / $r;

            $t = $pp / $d;


            if ($y < 0) {
                $this->sin = -$this->sin;
                $t = -$t;
            }

            $this->S[$k][$l] = 0;
            $this->update($k, -$t);
            $this->update($l, $t);

            /*
             * rotate rows and columns k and l
             */
            for ($i = 1; $i <= $k - 1; $i++) {
                $this->rotate($i, $k, $i, $l);
            }
            for ($i = $k + 1; $i <= $l - 1; $i++) {
                $this->rotate($k, $i, $i, $l);
            }
            for ($i = $l + 1; $i <= $this->n; $i++) {
                $this->rotate($k, $i, $l, $i);
            }

            /*
             * rotate eigenvectors
             */
            for ($i = 1; $i <= $this->n; $i++) {

                $rot_ik = $this->_E[$i][$k];
                $rot_il = $this->_E[$i][$l];

                $this->_E[$i][$k] = $this->cos * $rot_ik
                    - $this->sin * $rot_il;

                $this->_E[$i][$l] = $this->sin * $rot_ik
                    + $this->cos * $rot_il;
            }

            /*
             * rows k, l have changed, update rows ind_k, ind_l
             */
            $ind[$k] = $this->maxind($k);
            $ind[$l] = $this->maxind($l);

            $this->iter++;
        }
    }

    /**
     * Generate n x n identical matrix
     */
    protected final function I()
    {
        for ($i = 0; $i <= $this->n; $i++) {
            $this->_E[$i] = array_fill(0, $this->n + 1, 0);
            $this->_E[$i][$i] = 1;
        }
    }

    /**
     * Return index of largest off-diagonal element in row k
     *
     * @param $k
     * @return mixed
     */
    protected final function maxind($k)
    {
        $m = $k + 1;

        for ($i = $k + 2; $i <= $this->n; $i++) {
            if (abs($this->S[$k][$i]) > abs($this->S[$k][$m])) {
                $m = $i;
            }
        }

        return $m;
    }

    /**
     * Update e_k and its status
     */
    protected final function update($k, $t)
    {
        $y = $this->e[$k];
        $this->e[$k] = $y + $t;

        if ($this->changed[$k] && ($y == $this->e[$k])) {
            $this->changed[$k] = false;
            --$this->state;
        } else if (!$this->changed[$k] && ($y != $this->e[$k])) {
            $this->changed[$k] = true;
            ++$this->state;
        }
    }

    /**
     * Perform rotation of S_i,j and S_k,l
     */
    protected final function rotate($k, $l, $i, $j)
    {
        $rot_kl = $this->S[$k][$l];
        $rot_ij = $this->S[$i][$j];

        $this->S[$k][$l] = $this->cos * $rot_kl
            - $this->sin * $rot_ij;

        $this->S[$i][$j] = $this->sin * $rot_kl
            + $this->cos * $rot_ij;
    }

    /**
     * returned [V,e]
     * @return array
     */
    public final function eig()
    {
        $this->prepareEig();
        return $this->returning();
    }

    /**
     * Prepare V and e
     */
    protected final function prepareEig()
    {
        for ($i = 1; $i <= $this->n; $i++) {
            array_shift($this->_E[$i]);
        }
        array_shift($this->_E);
        array_shift($this->e);

        for ($i = 0; $i < $this->n; $i++) {
            for ($j = 0; $j < $this->n; $j++) {
                $this->V[$i][$j] = $this->_E[$j][$i];
            }
        }
    }

    /**
     * Return data
     * @return array
     */
    protected function returning()
    {
        return array($this->V, $this->e);
    }

    /**
     * Restore matrix S
     */
    public final function restore()
    {
        for ($k = 1; $k <= $this->n - 1; $k++) {
            for ($l = $k + 1; $l <= $this->n; $l++) {
                $this->S[$k][$l] = $this->S[$l][$k];
            }
        }
    }

    /**
     * Eigenvalues in descending order
     * @return array
     */
    public final function eigSort()
    {
        $this->prepareEig();
        $A = array();

        for ($i = 0; $i < $this->n; $i++) {
            $A[] = array('v' => $this->e[$i], 'p' => $i);
        }

        usort($A, array($this, 'sortFcn'));

        $e = array();
        $V = array();

        for ($i = 0; $i < $this->n; $i++) {
            $e[$i] = $this->e[$A[$i]['p']];
            $V[$i] = $this->V[$A[$i]['p']];
        }

        $this->e = $e;
        $this->V = $V;

        return $this->returning();
    }

    /**
     * Sort function for usort
     *
     * @param $a
     * @param $b
     * @return int
     */
    function sortFcn(&$a, &$b)
    {

        if ($a['v'] == $b['v'])
            return 0;
        return ($a['v'] > $b['v']) ? -1 : 1;
    }
}