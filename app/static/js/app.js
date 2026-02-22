(function () {
  const searchQueryInput = document.getElementById('search-query');
  const vcfSelect = document.getElementById('vcf-select');
  const searchBtn = document.getElementById('search-btn');
  const searchResultEl = document.getElementById('search-result');

  function escapeHtml(s) {
    const div = document.createElement('div');
    div.textContent = s == null ? '' : String(s);
    return div.innerHTML;
  }

  function setLoading(loading) {
    searchBtn.disabled = loading;
    if (loading) {
      searchResultEl.innerHTML = '<p class="flex items-center text-slate-500 text-sm"><span class="spinner" aria-hidden="true"></span> Searching…</p>';
    }
  }

  function setError(msg) {
    searchResultEl.innerHTML = '<div class="rounded-lg bg-red-50 border border-red-200 px-4 py-3 text-sm text-red-700">' + escapeHtml(msg) + '</div>';
  }

  function setMessage(text) {
    searchResultEl.innerHTML = '<pre class="rounded-lg bg-slate-50 border border-slate-200 p-4 text-sm text-slate-800 whitespace-pre-wrap font-mono leading-relaxed overflow-x-auto">' + escapeHtml(text) + '</pre>';
  }

  // Load VCF list on page load and init Select2
  function loadVcfList() {
    fetch('/api/vcf-list')
      .then(function (res) { return res.json(); })
      .then(function (data) {
        const list = data.vcf_list || [];
        vcfSelect.innerHTML = '<option value="">— Select a VCF —</option>' +
          list.map(function (item) {
            return '<option value="' + escapeHtml(item.path) + '">' + escapeHtml(item.label) + '</option>';
          }).join('');
        if (window.jQuery && window.jQuery.fn.select2) {
          window.jQuery('#vcf-select').select2({
            placeholder: '— Select a VCF —',
            allowClear: true,
            width: '100%'
          });
        }
      })
      .catch(function () {
        vcfSelect.innerHTML = '<option value="">— No VCFs configured (set VCF_DIR or VCF_PATHS) —</option>';
        if (window.jQuery && window.jQuery.fn.select2) {
          window.jQuery('#vcf-select').select2({
            placeholder: '— No VCFs configured —',
            allowClear: false,
            width: '100%'
          });
        }
      });
  }

  function runSearch() {
    const query = (searchQueryInput.value || '').trim();
    const vcfPath = (vcfSelect.value || '').trim();
    if (!query) {
      setError('Enter a gene symbol or location.');
      return;
    }
    if (!vcfPath) {
      setError('Select a VCF file.');
      return;
    }
    setLoading(true);
    fetch('/api/search', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: query, vcf_path: vcfPath })
    })
      .then(function (res) { return res.json().then(function (data) { return { res: res, data: data }; }); })
      .then(function (out) {
        setLoading(false);
        var res = out.res;
        var data = out.data;
        if (!res.ok) {
          var detail = data.detail || res.statusText;
          if (typeof detail === 'string') setError(detail);
          else setError(detail.msg || JSON.stringify(detail));
          return;
        }
        setMessage(data.message || '');
      })
      .catch(function (e) {
        setLoading(false);
        setError(e.message || 'Request failed.');
      });
  }

  searchBtn.addEventListener('click', runSearch);
  searchQueryInput.addEventListener('keydown', function (e) {
    if (e.key === 'Enter') runSearch();
  });

  loadVcfList();
})();
