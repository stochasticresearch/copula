#!/usr/bin/perl

use strict;
use warnings;
use Graph::Easy;

my $graph = Graph::Easy->new();

$graph->add_edge ('Z_{1}', 'X');
$graph->add_edge ('Z_{2}', 'X');
$graph->add_edge ('Z_{m}', 'X');
$graph->add_edge ('Z_{1}', 'Y');
$graph->add_edge ('Z_{2}', 'Y');
$graph->add_edge ('Z_{m}', 'Y');

$graph->add_node ('Z_{m+1}');
$graph->add_node ('Z_{n}');

print $graph->as_ascii( );              # prints:

print "\n\n\n";

my $graph2 = Graph::Easy->new();
$graph2->add_edge ('Z_{1}', 'X');
$graph2->add_edge ('Z_{2}', 'X');
$graph2->add_edge ('Z_{m}', 'X');
$graph2->add_edge ('Z_{1}', 'Y');
$graph2->add_edge ('Z_{2}', 'Y');
$graph2->add_edge ('Z_{m}', 'Y');
$graph2->add_edge ('E', 'X');
$graph2->add_edge ('E', 'Y');

$graph2->add_node ('Z_{m+1}');
$graph2->add_node ('Z_{n}');

print $graph2->as_ascii();

