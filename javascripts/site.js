
(function () {
  function applyLanguage(lang) {
    document.documentElement.lang = lang === 'zh' ? 'zh-CN' : 'en';
    document.querySelectorAll('[data-lang]').forEach(function (el) {
      el.style.display = el.getAttribute('data-lang') === lang ? '' : 'none';
    });
    document.querySelectorAll('[data-lang-switch]').forEach(function (btn) {
      btn.classList.toggle('active', btn.getAttribute('data-lang-switch') === lang);
    });
    localStorage.setItem('site-lang', lang);
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

  window.siteLang = {
    applyLanguage: applyLanguage,
    initLanguageToggle: initLanguageToggle
  };
})();
