#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);

my $root = shift @ARGV || '/root/autodl-tmp/faers-antithyroid-safety-signals';
my @default_quarters = qw(
  2014q4 2015q1 2015q2 2015q3 2015q4 2016q1 2016q2 2016q3 2016q4 2017q1
  2017q2 2017q3 2017q4 2018q1 2018q2 2018q3 2018q4 2019q1 2019q2 2019q3
  2019q4 2020q1 2020q2 2020q3 2020q4 2021q1 2021q2 2021q3 2021q4 2022q1
  2022q2 2022q3 2022q4 2023q1 2023q2 2023q3 2023q4 2024q1 2024q2 2024q3
  2024q4 2025q1
);
my @quarters = $ENV{FAERS_QUARTERS}
    ? grep { length } split(/[\s,]+/, $ENV{FAERS_QUARTERS})
    : @default_quarters;

my %primary_pt = (
    'ACUTE HEPATIC FAILURE' => 'HEPATIC',
    'HEPATIC FAILURE' => 'HEPATIC',
    'FULMINANT HEPATITIS' => 'HEPATIC',
    'TOXIC HEPATITIS' => 'HEPATIC',
    'DRUG-INDUCED LIVER INJURY' => 'HEPATIC',
    'LIVER INJURY' => 'HEPATIC',
    'HEPATOTOXICITY' => 'HEPATIC',
    'HEPATIC ENCEPHALOPATHY' => 'HEPATIC',
    'CHOLESTATIC HEPATITIS' => 'HEPATIC',
    'JAUNDICE' => 'HEPATIC',
    'AGRANULOCYTOSIS' => 'HEMATOLOGIC',
    'FEBRILE NEUTROPENIA' => 'HEMATOLOGIC',
    'NEUTROPENIA' => 'HEMATOLOGIC',
    'PANCYTOPENIA' => 'HEMATOLOGIC',
    'APLASTIC ANAEMIA' => 'HEMATOLOGIC',
    'APLASTIC ANEMIA' => 'HEMATOLOGIC',
    'BONE MARROW FAILURE' => 'HEMATOLOGIC',
    'MYELOSUPPRESSION' => 'HEMATOLOGIC',
);
my %serious_code = map { $_ => 1 } qw(DE LT HO DS CA OT);

sub split_faers {
    my ($line) = @_;
    chomp $line;
    $line =~ s/\r$//;
    return split(/\$/, $line, -1);
}

sub files_for {
    my ($q, $prefix) = @_;
    my $dir = "$root/raw/work/$q";
    opendir(my $dh, $dir) or die "cannot open $dir: $!";
    my @files;
    while (my $entry = readdir($dh)) {
        next if $entry eq '.' || $entry eq '..';
        my $path = "$dir/$entry";
        if (-d $path) {
            opendir(my $sub, $path) or next;
            while (my $child = readdir($sub)) {
                push @files, "$path/$child"
                    if $child =~ /^$prefix\d+Q[1-4](?:_new)?\.txt$/i;
            }
            closedir($sub);
        } elsif ($entry =~ /^$prefix\d+Q[1-4](?:_new)?\.txt$/i) {
            push @files, $path;
        }
    }
    closedir($dh);
    die "missing unique $prefix file for $q" unless @files == 1;
    return $files[0];
}

sub target_kind {
    my ($name, $ingredient) = @_;
    my $text = uc(($name // '') . ' ' . ($ingredient // ''));
    my $m = $text =~ /(?:^|[^A-Z])(?:METHIMAZOLE|THIAMAZOLE|TAPAZOLE|THYROZOL|THIAMAZOL|THACAPZOL|THYCAPZOL|TIRODRIL|STRUMAZOL|UNIMAZOLE|METHIZOL|METIMAZOL|METHEMAZOL|METIBASOL|THIMAZOLE)(?:[^A-Z]|$)/;
    my $p = $text =~ /(?:^|[^A-Z])(?:PROPYLTHIOURACIL|PROP[- ]?THIOURACIL|PROPHYLTHIOURACIL|PTU|PROPYCIL|TIURAGYL)(?:[^A-Z]|$)/;
    return 'MIXED_TARGETS' if $m && $p;
    return 'METHIMAZOLE' if $m;
    return 'PROPYLTHIOURACIL' if $p;
    return '';
}

sub date_days {
    my ($start, $event) = @_;
    return undef unless defined $start && defined $event;
    return undef unless $start =~ /^(\d{4})(\d{2})(\d{2})$/ && $event =~ /^(\d{4})(\d{2})(\d{2})$/;
    my ($sy, $sm, $sd) = ($1, $2, $3);
    my ($ey, $em, $ed) = ($event =~ /^(\d{4})(\d{2})(\d{2})$/);
    my $s = eval { require Time::Local; Time::Local::timelocal(0, 0, 0, $sd, $sm - 1, $sy) };
    my $e = eval { require Time::Local; Time::Local::timelocal(0, 0, 0, $ed, $em - 1, $ey) };
    return undef if !defined $s || !defined $e;
    return int(($e - $s) / 86400);
}

my (%latest, %case_count, %case_min_version, %case_max_version, %stats);
for my $q (@quarters) {
    my $file = files_for($q, 'DEMO');
    open(my $fh, '<', $file) or die "open $file: $!";
    <$fh>;
    while (my $line = <$fh>) {
        my @f = split_faers($line);
        next unless @f >= 8 && defined $f[1] && length $f[1];
        $stats{demo_rows}++;
        my ($primaryid, $caseid, $version) = @f[0, 1, 2];
        $version = 0 unless defined $version && $version =~ /^\d+$/;
        $case_count{$caseid}++;
        $case_min_version{$caseid} = $version if !defined $case_min_version{$caseid} || $version < $case_min_version{$caseid};
        $case_max_version{$caseid} = $version if !defined $case_max_version{$caseid} || $version > $case_max_version{$caseid};
        my $old = $latest{$caseid};
        if (!defined $old || $version > $old->{version} ||
            ($version == $old->{version} && $primaryid gt $old->{primaryid})) {
            $latest{$caseid} = {
                primaryid => $primaryid,
                version => $version,
                event_dt => $f[4] // '',
                age => $f[13] // '',
                age_grp => $f[15] // '',
                sex => $f[16] // '',
            };
        }
    }
    close $fh;
}

my %latest_pid = map { $latest{$_}{primaryid} => $_ } keys %latest;
my (%exp_ps, %exp_any, %target_seq, %drug_alias_counts);
for my $q (@quarters) {
    my $file = files_for($q, 'DRUG');
    open(my $fh, '<', $file) or die "open $file: $!";
    <$fh>;
    while (my $line = <$fh>) {
        my @f = split_faers($line);
        next unless @f >= 6;
        $stats{drug_rows}++;
        my ($primaryid, $seq, $role) = @f[0, 2, 3];
        next unless exists $latest_pid{$primaryid};
        my $kind = target_kind($f[4], $f[5]);
        next unless $kind;
        $drug_alias_counts{"${kind}_rows"}++;
        if (defined $role && ($role eq 'PS' || $role eq 'SS')) {
            $exp_any{$primaryid}{$kind} = 1;
            $target_seq{$primaryid}{$seq} = $kind if defined $seq;
        }
        if (defined $role && $role eq 'PS') {
            $exp_ps{$primaryid}{$kind} = 1;
        }
    }
    close $fh;
}

my (%reaction, %outcome, %therapy_start);
for my $q (@quarters) {
    my $file = files_for($q, 'REAC');
    open(my $fh, '<', $file) or die "open $file: $!";
    <$fh>;
    while (my $line = <$fh>) {
        my @f = split_faers($line);
        next unless @f >= 3 && exists $latest_pid{$f[0]};
        $stats{reac_rows}++;
        my $pt = uc($f[2] // '');
        $reaction{$f[0]}{$primary_pt{$pt}} = 1 if exists $primary_pt{$pt};
    }
    close $fh;

    my $outf = files_for($q, 'OUTC');
    open(my $ofh, '<', $outf) or die "open $outf: $!";
    <$ofh>;
    while (my $line = <$ofh>) {
        my @f = split_faers($line);
        next unless @f >= 3 && exists $latest_pid{$f[0]};
        $stats{outc_rows}++;
        $outcome{$f[0]}{$f[2]} = 1;
    }
    close $ofh;

    my $ther = files_for($q, 'THER');
    open(my $tfh, '<', $ther) or die "open $ther: $!";
    <$tfh>;
    while (my $line = <$tfh>) {
        my @f = split_faers($line);
        next unless @f >= 5 && exists $latest_pid{$f[0]};
        my ($pid, $seq, $start) = @f[0, 2, 3];
        next unless exists $target_seq{$pid}{$seq} && defined $start && length $start;
        $therapy_start{$pid} = $start if !defined $therapy_start{$pid} || $start lt $therapy_start{$pid};
    }
    close $tfh;
}

my %report;
for my $caseid (keys %latest) {
    my $pid = $latest{$caseid}{primaryid};
    my @ps = sort keys %{ $exp_ps{$pid} // {} };
    my @any = sort keys %{ $exp_any{$pid} // {} };
    my $ps_group = @ps == 1 ? $ps[0] : (@ps > 1 ? 'MIXED_TARGETS' : 'NEITHER');
    my $any_group = @any == 1 ? $any[0] : (@any > 1 ? 'MIXED_TARGETS' : 'NEITHER');
    my $serious = 0;
    $serious = 1 if grep { exists $outcome{$pid}{$_} } keys %serious_code;
    my $days = date_days($therapy_start{$pid}, $latest{$caseid}{event_dt});
    $report{$caseid} = {
        pid => $pid,
        ps_group => $ps_group,
        any_group => $any_group,
        hepatic => $reaction{$pid}{HEPATIC} ? 1 : 0,
        hematologic => $reaction{$pid}{HEMATOLOGIC} ? 1 : 0,
        any_target_outcome => ($reaction{$pid}{HEPATIC} || $reaction{$pid}{HEMATOLOGIC}) ? 1 : 0,
        serious => $serious,
        sex => $latest{$caseid}{sex} // '',
        age_grp => $latest{$caseid}{age_grp} // '',
        onset_days => $days,
    };
}

sub ror_line {
    my ($outcome_key, $group_key, $target, $background, $stratum, $filter_key, $filter_value) = @_;
    my ($a, $b, $c, $d) = (0, 0, 0, 0);
    for my $r (values %report) {
        next if defined $filter_key && $r->{$filter_key} ne $filter_value;
        my $group = $r->{$group_key};
        next unless $group eq $target || $group eq $background;
        my $event = $r->{$outcome_key} ? 1 : 0;
        if ($group eq $target) { $event ? $a++ : $b++; }
        else { $event ? $c++ : $d++; }
    }
    my ($aa, $bb, $cc, $dd) = map { $_ + 0.5 } ($a, $b, $c, $d);
    my $ror = ($aa * $dd) / ($bb * $cc);
    my $se = sqrt(1 / $aa + 1 / $bb + 1 / $cc + 1 / $dd);
    my ($lo, $hi) = (exp(log($ror) - 1.96 * $se), exp(log($ror) + 1.96 * $se));
    return join("\t", $stratum, $outcome_key, $group_key, $target, $background,
        $a, $b, $c, $d, sprintf('%.6f', $ror), sprintf('%.6f', $lo), sprintf('%.6f', $hi));
}

make_path("$root/results");
open(my $ror, '>', "$root/results/full_window_ror.tsv") or die $!;
print $ror "stratum\toutcome\texposure_definition\ttarget\tbackground\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n";
my @strata = (['ALL', undef, undef], ['SEX_F', 'sex', 'F'], ['SEX_M', 'sex', 'M']);
for my $age (qw(A E N P T)) { push @strata, ["AGE_$age", 'age_grp', $age]; }
for my $s (@strata) {
    for my $outcome_key (qw(hepatic hematologic any_target_outcome serious)) {
        for my $target (qw(METHIMAZOLE PROPYLTHIOURACIL)) {
            print $ror ror_line($outcome_key, 'ps_group', $target, 'NEITHER', @$s), "\n";
        }
    }
}
close $ror;

sub group_counts {
    my ($key) = @_;
    my %counts;
    $counts{$report{$_}{$key}}++ for keys %report;
    return %counts;
}

open(my $sum, '>', "$root/results/full_window_summary.txt") or die $!;
print $sum "FAERS antithyroid analysis\nwindow=" . join(',', @quarters) . "\nquarters=" . scalar(@quarters) . "\nreport_unit=latest caseversion per caseid; primaryid tiebreak\nprimary_exposure=role_cod PS\nsensitivity_exposure=role_cod PS or SS\nbackground=NEITHER target drug\noutcome=explicit frozen PT list\nror=0.5 continuity correction; Wald 95% CI\n\nraw_rows\n";
print $sum join("\n", map { "$_=$stats{$_}" } sort keys %stats), "\n\ncase_deduplication\n";
my $repeated = grep { $case_count{$_} > 1 } keys %case_count;
my $extra = 0; $extra += $case_count{$_} - 1 for keys %case_count;
my $versioned = grep { $case_min_version{$_} < $case_max_version{$_} } keys %case_count;
print $sum "latest_caseids=" . scalar(keys %latest) . "\nrepeated_caseids=$repeated\nextra_rows_removed=$extra\ncaseids_with_version_change=$versioned\n\nps_groups\n";
my %ps_counts = group_counts('ps_group');
print $sum join("\n", map { "$_=$ps_counts{$_}" } sort keys %ps_counts), "\n\nps_or_ss_groups\n";
my %any_counts = group_counts('any_group');
print $sum join("\n", map { "$_=$any_counts{$_}" } sort keys %any_counts), "\n\ntarget_alias_rows\n";
print $sum join("\n", map { "$_=$drug_alias_counts{$_}" } sort keys %drug_alias_counts), "\n\nvalid_onset_days=" . scalar(grep { defined $report{$_}{onset_days} && $report{$_}{onset_days} >= 0 } keys %report) . "\nror_file=results/full_window_ror.tsv\n";
close $sum;
print "WROTE full_window_summary.txt\nWROTE full_window_ror.tsv\n";
