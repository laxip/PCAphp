# PCAphp
Principal Component Analysis PHP

##Example

```
include( "pca.php" );

$points = [
    [2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2, 1, 1.5, 1.1],
    [2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9]
];

$p = new PCA\PCA($points);
$p->changeDimension(1);
$p->applayingPca();
print_r($p->getNewData());
print_r($p->getSimilarity());
```
