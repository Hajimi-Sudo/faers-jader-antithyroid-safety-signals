#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);

my $root = $ENV{FAERS_PROJECT_ROOT} || '/root/autodl-tmp/faers-antithyroid-safety-signals';
my $case = "$root/literature_case_evidence";
my $out = "$root/submission_qc";
make_path($out);

sub count_data_rows {
    my ($path) = @_;
    open my $fh, '<', $path or die "cannot read $path: $!";
    my $n = -1;
    $n++ while <$fh>;
    close $fh;
    return $n;
}

sub count_matching_rows {
    my ($path, $pattern) = @_;
    open my $fh, '<', $path or die "cannot read $path: $!";
    my $n = 0;
    while (my $line = <$fh>) { $n++ if $line =~ $pattern; }
    close $fh;
    return $n;
}

sub field_index {
    my ($path, $name) = @_;
    open my $fh, '<', $path or die "cannot read $path: $!";
    my $header = <$fh>;
    close $fh;
    chomp $header;
    my @fields = split /\t/, $header, -1;
    for my $i (0 .. $#fields) { return $i if $fields[$i] eq $name; }
    die "missing field $name in $path";
}

sub count_field_value {
    my ($path, $field, $wanted) = @_;
    my $index = field_index($path, $field);
    open my $fh, '<', $path or die "cannot read $path: $!";
    <$fh>;
    my $n = 0;
    while (my $line = <$fh>) {
        chomp $line;
        my @values = split /\t/, $line, -1;
        $n++ if defined $values[$index] && $values[$index] eq $wanted;
    }
    close $fh;
    return $n;
}

my $screen_rows = count_data_rows("$case/screening_audit.tsv");
my $curated_rows = count_data_rows("$case/curated_case_literature_summary.tsv");
my $manual_rows = count_data_rows("$case/manual_case_review_sheet.tsv");
my $included = count_field_value("$case/screening_audit.tsv", 'include_for_manual_case_review', 'true');
my $direct = count_field_value("$case/screening_audit.tsv", 'evidence_tier', 'direct_case_signal');
my $pending = count_matching_rows("$case/manual_case_review_sheet.tsv", qr/\tmanual_verification_pending\t/);

open my $flow, '>', "$out/case_evidence_flow.tsv" or die $!;
print $flow "stage\tn\tinterpretation\n";
print $flow "PubMed/Europe PMC candidate articles\t$screen_rows\tsearch and metadata retrieval\n";
print $flow "Manual case-review candidates\t$included\theuristic screen; not final inclusion\n";
print $flow "Direct-attribution screen candidates\t$direct\tautomated attribution signal\n";
print $flow "Curated publication-level evidence units\t$curated_rows\tcontextual summary\n";
print $flow "Author-verified evidence units\t" . ($manual_rows - $pending) . "\trequired before submission\n";
close $flow;

open my $map, '>', "$out/READUS_PV_MAPPING.md" or die $!;
print $map <<'EOF';
# READUS-PV Submission Mapping

Status: BLOCKED pending author-level clinical verification and completion of the author declarations.

| READUS-PV domain | Evidence artifact | Status | Remaining action |
|---|---|---|---|
| Background and objective | manuscript.tex, Introduction | present | retain reporting-signal wording |
| Study design and data source | manuscript.tex, Methods 2.1 | present | add final access dates/release identifiers |
| Study population and report unit | final_qc_report.json; full_window_summary.txt | present | retain database-specific duplicate rules |
| Exposure definition | analysis scripts; manuscript.tex, Methods 2.3 | present | place full alias dictionary in supplement |
| Outcome and serious-outcome definitions | manuscript.tex, Methods 2.4; OT audit outputs | present after correction | keep DE/LT/HO/DS/CA/OT and explicitly exclude RI |
| Statistical analysis | aggregate ROR/PRR/IC outputs; Methods 2.5 | present | add exact software/version details |
| Results flow and exclusions | case_evidence_flow.tsv; population table | partial | reconcile final author-verified case flow |
| Disproportionality results | supplement_ror_prr.tsv; supplement_ic.tsv | present | retain denominators and uncertainty intervals |
| Case-by-case clinical evidence | manual_case_review_sheet.tsv; curated summary | partial, critical | verify each included publication and mark inclusion/exclusion |
| Interpretation and limitations | manuscript.tex, Discussion | present | add publication-selection bias after verification |
| Data/code/declarations | manuscript.tex, Declarations | partial | complete funding, COI, CRediT, protocol and code statements |

The case-evidence screen is a review aid. Its counts are publication-level
evidence units and are not independent patient counts, a systematic review, a
causality adjudication, or validation of FAERS/JADER.
EOF
close $map;

open my $qc, '>', "$out/submission_qc.json" or die $!;
print $qc "{\n";
print $qc "  \"status\": \"blocked_pending_author_case_verification\",\n";
print $qc "  \"raw_data_read\": false,\n  \"report_level_matrix_read\": false,\n  \"faers_final_qc_present\": true,\n";
print $qc "  \"case_screen\": {\"candidate_articles\": $screen_rows, \"manual_case_review_candidates\": $included, \"direct_attribution_screen_candidates\": $direct, \"curated_publication_units\": $curated_rows, \"manual_review_rows\": $manual_rows, \"manual_verification_pending\": $pending},\n";
print $qc "  \"blocking_items\": [\"author-level verification of each included publication-level evidence unit\", \"completion of READUS-PV checklist sign-off\", \"completion of author declarations\", \"final reference identity and DOI verification\"],\n";
print $qc "  \"claim_boundary\": \"reporting signals and contextual published-case evidence only; no incidence, causality, pooled risk, or confirmed interaction\"\n}\n";
close $qc;
print "wrote submission QC: candidates=$screen_rows curated=$curated_rows pending_manual=$pending\n";
