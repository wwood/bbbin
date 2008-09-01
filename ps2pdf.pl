#!/usr/bin/perl -w

`ls *.ps>/tmp/ps2pdftmp`;
open FILES, "/tmp/ps2pdftmp";
foreach (<FILES>){
  chomp;
  `ps2pdf $_`;
}
