<?php namespace PCA;

include( "Math/jacobiEigenvalueAlgorithm.php" );

/**
 * Class for calculate Principal Component Analysis
 *
 * @author laxip
 *
 * Class PCA
 * @package PCA
 */
class PCA
{
    /**
     * length of column for displayed matrix
     */
    const PRINT_LENGTH = 12;

    /**
     * @var array data matrix
     */
    protected $A;

    /**
     * @var array mean vector
     */
    protected $M;

    /**
     * @var array normalized data matrix
     */
    protected $A_;

    /**
     * @var array covariance matrix
     */
    protected $C;

    /**
     * @var array eigenvectors
     */
    protected $Z_;

    /**
     * @var array data in new space
     */
    protected $X_;

    /**
     * @var array eigenvalues
     */
    protected $L;

    /**
     * @var array cropped data in nominal space
     */
    protected $newA;

    /**
     * @var null similarity
     */
    protected $s = null;

    /**
     * @var int new space dimension
     */
    protected $p;

    /**
     * PCA constructor.
     * @param $_A
     * @param null $_n
     * @param null $_m
     */
    function __construct($_A, $_n = null, $_m = null)
    {
        $this->A = $_A;

        $this->m = (is_null($_m)) ? count($this->A) : $_m; //number of columns
        $this->n = (is_null($_n)) ? count($this->A[0]) : $_n; //number of rows

        $this->prepareCovMatrix();
    }

    /**
     * Make steps
     */
    protected final function prepareCovMatrix()
    {
        $this->mean();
        $this->makeA_();
        $this->makeCov();

        $eig = new Math\jacobiEigenvalueAlgorithm($this->C, $this->m);
        list($this->Z_, $this->L) = $eig->eigSort();
    }

    /**
     * Function calculate mean value for each column
     */
    protected final function mean()
    {
        for ($i = 0; $i < $this->m; $i++) {
            $sum = 0;

            for ($j = 0; $j < $this->n; $j++) {
                $sum += $this->A[$i][$j];
            }
            $this->M[$i] = $sum / $this->n;
        }
    }

    /**
     * Calculate normalized matrix. For each column subtract mean value
     */
    protected final function makeA_()
    {
        for ($i = 0; $i < $this->m; $i++) {
            for ($j = 0; $j < $this->n; $j++) {
                $this->A_[$i][$j] = $this->A[$i][$j] - $this->M[$i];
            }
        }
    }

    /**
     * Calcaulate covariance matrix
     */
    protected final function makeCov()
    {
        for ($i = 0; $i < $this->m; $i++) {
            for ($j = 0; $j <= $i; $j++) {
                $c = $this->cov($i, $j);
                $this->C[$i][$j] = $this->C[$j][$i] = $c;
            }
        }
    }

    /**
     * Calculate covariance between two columns
     */
    protected final function cov($i, $j)
    {
        $c = 0;
        for ($k = 0; $k < $this->n; $k++) {
            $c += $this->A_[$i][$k] * $this->A_[$j][$k];
        }

        return $c / ($this->n - 1);
    }

    /**
     * Print matrix by column index
     */
    public function printViaColumn(array &$a)
    {
        $m = count($a);
        $n = count($a[0]);

        for ($j = 0; $j < $n; $j++) { //row
            for ($i = 0; $i < $m; $i++) { //column

                $sh = substr($a[$i][$j], 0, self::PRINT_LENGTH);
                $sh = str_pad($sh, self::PRINT_LENGTH, " ", STR_PAD_LEFT);

                echo $sh . (($i < $m - 1) ? " | " : "");
            }
            echo PHP_EOL;
        }
    }

    /**
     * Print matrix by row index
     */
    public function printViaRow(array &$a)
    {

    }

    /**
     * Set new dimension for space
     */
    public function changeDimension($_p)
    {
        $this->p = min($this->m, $_p);
        $this->Z_ = array_slice($this->Z_, 0, $this->p);

        $this->makeX_();
    }

    /**
     * Calculate data in new space
     */
    protected final function makeX_()
    {
        for ($i = 0; $i < $this->p; $i++) {
            for ($j = 0; $j < $this->n; $j++) {
                $sum = 0;

                for ($k = 0; $k < $this->m; $k++) {
                    $sum += $this->A_[$k][$j] * $this->Z_[$i][$k];
                }

                $this->X_[$i][$j] = $sum;
            }
        }
    }

    /**
     * Make obtained data in nominal space
     */
    public function applayingPca()
    {
        for ($i = 0; $i < $this->m; $i++) {
            for ($j = 0; $j < $this->n; $j++) {
                $sum = 0;

                for ($k = 0; $k < $this->p; $k++) {
                    $sum += $this->X_[$k][$j] * $this->Z_[$k][$i];
                }

                $this->newA[$i][$j] = $sum;
            }
        }

        $this->restore();
    }

    /**
     * Add mean for each element
     */
    protected final function restore()
    {
        for ($i = 0; $i < $this->m; $i++) {
            for ($j = 0; $j < $this->n; $j++) {
                $this->newA[$i][$j] += $this->M[$i];
            }
        }
    }

    /**
     * Return Principal Components
     * @return array
     */
    public function getPC()
    {
        return $this->X_;
    }

    /**
     * Return projected data
     * @return array
     */
    public function getProjectedData()
    {
        return $this->X_;
    }

    /**
     * Return obtained data in nominal space
     * @return array
     */
    public function getNewData()
    {
        return $this->newA;
    }

    /**
     * Return similarity
     * @return float
     */
    function getSimilarity()
    {
        if (!is_null($this->s))
            return $this->s;

        $upp = 0;
        for ($i = 0; $i < $this->p; $i++) {
            $upp += $this->L[$i];
        }

        $upm = $upp;
        for (; $i < $this->m; $i++) {
            $upm += $this->L[$i];
        }

        $this->s = $upp / $upm;
        return $this->s;
    }
}
