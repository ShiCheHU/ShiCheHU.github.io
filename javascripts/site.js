/**
 * Shicheng Hu — site.js
 * Language toggle + blog-post markdown renderer
 * Zero dependencies, ~1 KB
 */
(function () {
  /* ========== Language Toggle ========== */
  function applyLanguage(lang) {
    document.documentElement.lang = lang === 'zh' ? 'zh-CN' : 'en';
    document.querySelectorAll('[data-lang]').forEach(function (el) {
      el.style.display = el.getAttribute('data-lang') === lang ? '' : 'none';
    });
    document.querySelectorAll('[data-lang-switch]').forEach(function (btn) {
      btn.classList.toggle('active', btn.getAttribute('data-lang-switch') === lang);
    });
    localStorage.setItem('site-lang', lang);
    // Update lang-aware links
    document.querySelectorAll('a[data-lang-link-en][data-lang-link-zh]').forEach(function (a) {
      a.setAttribute('href', lang === 'zh' ? a.getAttribute('data-lang-link-zh') : a.getAttribute('data-lang-link-en'));
    });
  }

  function initLanguageToggle(defaultLang) {
    var lang = localStorage.getItem('site-lang') || defaultLang || 'en';
    applyLanguage(lang);
    document.querySelectorAll('[data-lang-switch]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        applyLanguage(btn.getAttribute('data-lang-switch'));
      });
    });
  }

  /* ========== Blog Post Renderer ========== */
  function initBlogPost() {
    var container = document.getElementById('post-content');
    if (!container) return;

    var params = new URLSearchParams(window.location.search);
    var lang = params.get('lang') || localStorage.getItem('site-lang') || 'zh';
    var post = params.get('post');
    var file = params.get('file');

    var url;
    if (post) {
      // Standard blog: blog/posts/{lang}/{post}.md
      url = 'blog/posts/' + lang + '/' + post + '.md';
    } else if (file) {
      // Tech note: append .en.md for English
      url = lang === 'en' ? file.replace(/\.md$/, '.en.md') : file;
    } else {
      container.innerHTML = '<p style="color:var(--c-muted)">No post specified.</p>';
      return;
    }

    fetch(url)
      .then(function (res) {
        if (!res.ok) throw new Error('Not found: ' + url);
        // Force UTF-8 decoding to prevent garbled Chinese on some servers
        return res.arrayBuffer();
      })
      .then(function (buffer) {
        return new TextDecoder('utf-8').decode(buffer);
      })
      .then(function (md) {
        container.innerHTML = marked.parse(md);
        document.title = (container.querySelector('h1') || {}).textContent || 'Blog';
        // KaTeX rendering
        if (typeof renderMathInElement === 'function') {
          renderMathInElement(container, {
            delimiters: [
              { left: '$$', right: '$$', display: true },
              { left: '$', right: '$', display: false },
              { left: '\\[', right: '\\]', display: true },
              { left: '\\(', right: '\\)', display: false }
            ],
            throwOnError: false
          });
        }
      })
      .catch(function (err) {
        container.innerHTML = '<p style="color:#cf222e">Error: ' + err.message + '</p>';
      });
  }

  /* ========== Language switch with page reload (for blog-post) ========== */
  function initBlogLangSwitch() {
    document.querySelectorAll('[data-lang-switch]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var nextLang = btn.getAttribute('data-lang-switch');
        localStorage.setItem('site-lang', nextLang);
        // If on blog-post page, reload with new lang param
        var isBlogPost = document.getElementById('post-content');
        if (isBlogPost) {
          var u = new URL(window.location.href);
          u.searchParams.set('lang', nextLang);
          window.location.href = u.toString();
        }
        // For tech notes (?file=), just re-fetch
        var u = new URL(window.location.href);
        if (u.searchParams.has('file') && !u.searchParams.has('post')) {
          u.searchParams.set('lang', nextLang);
          window.location.href = u.toString();
        }
      });
    });
  }

  /* ========== Accordion (personal page) ========== */
  function initAccordion() {
    document.querySelectorAll('.accordion-btn').forEach(function (btn) {
      btn.addEventListener('click', function () {
        this.classList.toggle('open');
        var panel = this.nextElementSibling;
        if (panel) panel.classList.toggle('show');
      });
    });
  }

  /* ========== Exports ========== */
  window.siteLang = {
    applyLanguage: applyLanguage,
    initLanguageToggle: initLanguageToggle
  };

  /* ========== Auto-init ========== */
  document.addEventListener('DOMContentLoaded', function () {
    initBlogPost();
    initBlogLangSwitch();
    initAccordion();
  });
})();
