#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);

my $root = $ENV{FAERS_PROJECT_ROOT} || '/root/autodl-tmp/faers-antithyroid-safety-signals';
my $base = "$root/literature_case_evidence";
my $outdir = "$root/submission_qc";
make_path($outdir);

sub clean {
    my ($x) = @_;
    $x = '' unless defined $x;
    $x =~ s/<[^>]+>/ /g;
    $x =~ s/&amp;/\&/g; $x =~ s/&lt;/</g; $x =~ s/&gt;/>/g;
    $x =~ s/&#x([0-9A-Fa-f]+);/chr(hex($1))/ge;
    $x =~ s/&#([0-9]+);/chr($1)/ge;
    $x =~ s/\s+/ /g;
    $x =~ s/^\s+|\s+$//g;
    return $x;
}

sub norm {
    my ($x) = @_;
    $x = lc clean($x);
    $x =~ s/[^a-z0-9]+//g;
    return $x;
}

sub read_tsv {
    my ($path) = @_;
    open my $fh, '<', $path or die "cannot read $path: $!";
    my $header = <$fh>; chomp $header;
    my @fields = split /\t/, $header, -1;
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        my @values = split /\t/, $line, -1;
        my %row; @row{@fields} = @values;
        push @rows, \%row;
    }
    close $fh;
    return \@rows;
}

sub strip_xml_body {
    my ($path) = @_;
    return '' unless -f $path;
    open my $fh, '<', $path or return '';
    local $/; my $xml = <$fh>; close $fh;
    my ($body) = $xml =~ m{<body\b[^>]*>(.*?)</body>}is;
    return clean($body || $xml);
}

sub first_match {
    my ($re, $text) = @_;
    return $1 if $text =~ /$re/is;
    return '';
}

my @curated = @{read_tsv("$base/curated_case_literature_summary.tsv")};
my @candidates = @{read_tsv("$base/pubmed_candidates.tsv")};
my %candidate = map { $_->{pmid} => $_ } @candidates;

my @pmids = map { $_->{pmid} } @curated;
my $idlist = join(',', @pmids);
my $efetch = "/tmp/faers_pubmed_verify_$$.xml";
my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$idlist&retmode=xml";
system("curl -L --retry 2 --connect-timeout 20 --max-time 120 -sS '$url' -o '$efetch'") == 0 or die "PubMed efetch failed\n";
open my $xfh, '<', $efetch or die $!; local $/; my $xml = <$xfh>; close $xfh; unlink $efetch;

my %pub;
while ($xml =~ m{<PubmedArticle>(.*?)</PubmedArticle>}gis) {
    my $article = $1;
    my $pmid = first_match(qr{<PMID[^>]*>(.*?)</PMID>}, $article);
    next unless $pmid;
    my $title = first_match(qr{<ArticleTitle[^>]*>(.*?)</ArticleTitle>}, $article);
    my $year = first_match(qr{<PubDate>.*?<Year>(\d{4})</Year>}s, $article);
    $year ||= first_match(qr{<PubDate>.*?<MedlineDate>.*?(\d{4})}s, $article);
    my $doi = '';
    while ($article =~ m{<ArticleId IdType="doi">(.*?)</ArticleId>}gis) { $doi = clean($1); last; }
    my @abstract;
    while ($article =~ m{<AbstractText[^>]*>(.*?)</AbstractText>}gis) { push @abstract, clean($1); }
    $pub{$pmid} = {title => clean($title), year => $year, doi => $doi, abstract => join(' ', @abstract)};
}

my @verified;
for my $row (@curated) {
    my $p = $pub{$row->{pmid}} || {};
    my $c = $candidate{$row->{pmid}} || {};
    my $title_match = norm($p->{title}) eq norm($row->{title}) ? 'verified' : 'mismatch_or_unavailable';
    my $year_match = ($p->{year} && $row->{year} && $p->{year} eq $row->{year}) ? 'verified' : 'mismatch_or_unavailable';
    my $doi_status = 'not_available_in_pubmed';
    if ($row->{doi} && $p->{doi} && lc($row->{doi}) eq lc($p->{doi})) { $doi_status = 'verified'; }
    elsif ($row->{doi} && $p->{doi}) { $doi_status = 'mismatch'; }
    elsif ($row->{doi} && !$p->{doi}) { $doi_status = 'curated_doi_not_returned_by_efetch'; }

    my $full = strip_xml_body("$base/fulltext_xml/$row->{pmcid}.xml");
    my $source = clean(join(' ', $p->{title} || '', $p->{abstract} || '', $full));
    my $drug_re = $row->{target_drug} eq 'methimazole' ? qr{\b(?:methimazole|thiamazole)\b}i : qr{\b(?:propylthiouracil|PTU)\b}i;
    my $outcome_re = $row->{outcome_group} eq 'hepatic' ? qr{hepat|liver|jaundice|cholest|bilirubin|transaminase}i : qr{agranul|neutrop|pancyt|aplastic|leukopen|marrow|granulocyt}i;
    my $drug_event = ($source =~ /$drug_re/ && $source =~ /$outcome_re/) ? 'verified_in_public_text' : 'not_confirmed';
    my $course = ($source =~ /(?:case|patient|woman|man|girl|boy|child)/i && $source =~ /(?:discontinu|withdraw|stop|recover|improv|resolv|death|died|fatal|transplant|hospital|rechalleng)/i) ? 'verified_in_public_text' : 'unclear';
    my $alternative = ($source =~ /(?:viral|autoimmune|infection|sepsis|biliary|obstruction|thyroid|pregnan|alcohol|concomitant|other drug|hepatitis [abc])/i) ? 'flagged_or_requires_review' : 'not_identified_in_public_text';
    my $source_level = $full ? 'PubMed abstract plus accessible Europe PMC full text' : 'PubMed abstract; full text unavailable in remote corpus';
    my $evidence_status = ($title_match eq 'verified' && $year_match eq 'verified' && $drug_event eq 'verified_in_public_text' && $course eq 'verified_in_public_text') ? ($full ? 'public_text_verified_fulltext' : 'public_text_verified_abstract_only') : 'identity_or_clinical_evidence_incomplete';
    push @verified, {
        pmid => $row->{pmid}, pmcid => $row->{pmcid}, curated_doi => $row->{doi}, verified_doi => $p->{doi} || '',
        curated_title => $row->{title}, pubmed_title => $p->{title} || '', year => $row->{year}, pubmed_year => $p->{year} || '',
        title_identity => $title_match, year_identity => $year_match, doi_identity => $doi_status,
        target_drug => $row->{target_drug}, outcome_group => $row->{outcome_group}, drug_event_evidence => $drug_event,
        clinical_course_evidence => $course, competing_cause_status => $alternative, source_level => $source_level,
        evidence_status => $evidence_status, author_signoff => 'pending',
        note => $full ? 'Identity and clinical evidence checked against PubMed metadata/abstract and accessible full text.' : 'Identity and clinical evidence checked against PubMed metadata/abstract; full-text-specific fields remain limited.'
    };
}

my @fields = qw(pmid pmcid curated_doi verified_doi curated_title pubmed_title year pubmed_year title_identity year_identity doi_identity target_drug outcome_group drug_event_evidence clinical_course_evidence competing_cause_status source_level evidence_status author_signoff note);
open my $out, '>', "$outdir/case_evidence_verification.tsv" or die $!;
print $out join("\t", @fields), "\n";
for my $r (@verified) { print $out join("\t", map { ($r->{$_} // '') =~ s/[\t\r\n]+/ /gr } @fields), "\n"; }
close $out;

my $full_n = scalar grep { $_->{evidence_status} eq 'public_text_verified_fulltext' } @verified;
my $abs_n = scalar grep { $_->{evidence_status} eq 'public_text_verified_abstract_only' } @verified;
my $incomplete_n = scalar grep { $_->{evidence_status} eq 'identity_or_clinical_evidence_incomplete' } @verified;
open my $sum, '>', "$outdir/case_evidence_verification_summary.json" or die $!;
print $sum "{\n  \"status\": \"public_text_verification_complete\",\n  \"records\": " . scalar(@verified) . ",\n  \"public_text_verified_fulltext\": $full_n,\n  \"public_text_verified_abstract_only\": $abs_n,\n  \"identity_or_clinical_evidence_incomplete\": $incomplete_n,\n  \"author_signoff_pending\": " . scalar(@verified) . ",\n  \"claim_boundary\": \"publication-level public-text verification; not independent patient counts, causal adjudication, or systematic-review completeness\"\n}\n";
close $sum;
print "verified_records=" . scalar(@verified) . " fulltext=$full_n abstract_only=$abs_n incomplete=$incomplete_n\n";
