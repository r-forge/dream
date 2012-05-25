
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

For information on how to use dream, please run in R:<ul>
<li><i>help("dreamCalibrate")</i> - to calibrate a function using dream</li>
<li><i>help("dream")</i> - for low-level interface </li>
<li><i>demo(example1)</i> - Fitting a banana shaped distribution</li>
<li><i>demo(example2)</i> - Fitting an n-dimensional Gaussian distribution</li>
<li><i>demo(FME.nonlinear.model)</i> - Calibrating the non-linear model shown
      in the <a href="http://fme.r-forge.r-project.org/">FME package</a> vignette</li>
<li><i>demo(FME.nonlinear.model_parallelisation)</i> - Example of parallelisation using the <a href="http://cran.r-project.org/web/packages/snow/index.html">SNOW package</a></li>
<li><i>demo(parallelisation_chain_id)</i> - Example of parallelisation when DREAM calls an external model using batch files in separate folders</li>
</ul>
<p>
To cite the DREAM algorithm please use:
<ul>
<li>
  Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Robinson, B. A.,
  Hyman, J. M., Higdon, D., 2009. Accelerating Markov chain Monte Carlo
  simulation by differential evolution with self-adaptive randomized
  subspace sampling. <b>International Journal of Nonlinear Sciences
  and Numerical Simulation</b> 10 (3), 273-290. DOI: <a href="http://dx.doi.org/10.1515/IJNSNS.2009.10.3.273">10.1515/IJNSNS.2009.10.3.273</a>
</li></ul>
<p>
To cite the dream package, please use:<br><ul><li>
 Joseph Guillaume and Felix Andrews (2012). dream: DiffeRential
  Evolution Adaptive Metropolis. R package version 0.4-2. URL
  <a href="http://CRAN.R-project.org/package=dream">http://CRAN.R-project.org/package=dream</a>
</li></ul>

For additional information on the algorithm also see:<ul>
<li>
  Vrugt, J. A., ter Braak, C. J. F., Gupta, H. V., Robinson, B. A.,
  2009. Equifinality of formal (DREAM) and informal (GLUE) Bayesian
  approaches in hydrologic modeling?
  <b>Stochastic Environmental Research and Risk Assessment</b> 23 (7), 1011--1026.
  DOI: <a href="http://dx.doi.org/10.1007/s00477-008-0274-y">10.1007/s00477-008-0274-y</a>
</li>
</ul>

<p>
This implementation of DREAM has been tested against the original Matlab implementation. See <a href="./matlab_test/example1.R">example1.R</a> and <a href="./matlab_test/example2.R">example2.R</a>
<p>
Please note that the dream_zs and dream_d algorithms may be superior in your circumstances. These are not implemented in this package. Please read the following references for details:
<ul>
<li>
Vrugt, J. A. and Ter Braak, C. J. F. (2011) DREAM(D): an adaptive Markov Chain Monte Carlo simulation algorithm to solve discrete, noncontinuous, and combinatorial posterior parameter estimation problems, <b>Hydrol. Earth Syst. Sci.</b>, 15, 3701-3713, DOI: <a href="http://dx.doi.org/10.5194/hess-15-3701-2011">10.5194/hess-15-3701-2011</a>
</li>
<li>
ter Braak, C. and J. Vrugt (2008). Differential Evolution Markov Chain
with snooker updater and fewer chains. <b>Statistics and Computing</b> 18(4): 435-446 DOI: <a href="http://dx.doi.org/10.1007/s11222-008-9104-9">10.1007/s11222-008-9104-9</a>
</li>
<li>
Laloy,E., and J.A. Vrugt.  2012. High-dimensional posterior exploration
of hydrologic models using multiple-try DREAM(ZS) and high-performance
computing.  <b>Water Resources Research</b>, 48, W0156. DOI <a href="http://dx.doi.org/10.1029/2011WR010608">10.1029/2011WR010608</a>
</li></ul>

</body>
</html>
