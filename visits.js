/* ===== Visit log (shared by both labs) =====
   One record per visit (once per browser tab session) in Firestore visits/: time, site, page, device type,
   operating system, browser, screen, language, time zone, referrer, IP, and the account if the visitor is signed in.
   Only admins can read it (see firestore.rules). Runs only on the public GitHub Pages copy. ===== */
(function () {
  'use strict';
  if (!/\.github\.io$/.test(location.hostname)) return;
  var APP = /psych-stats-lab/.test(location.pathname) ? 'stats' : 'scales';
  var KEY = 'visit-logged-' + APP;
  try { if (sessionStorage.getItem(KEY)) return; sessionStorage.setItem(KEY, '1'); } catch (e) { }
  var cut = function (v, n) { return String(v == null ? '' : v).slice(0, n); };

  function device(ua) {
    var tablet = /iPad|Tablet|PlayBook|Silk|(Android(?!.*Mobile))/i.test(ua) || (/Macintosh/.test(ua) && navigator.maxTouchPoints > 1);
    var mobile = !tablet && /Mobi|iPhone|iPod|Android.*Mobile|Windows Phone/i.test(ua);
    var dev = tablet ? 'تابلت' : mobile ? 'موبايل' : 'كمبيوتر';
    var m, os = 'غير معروف';
    if ((m = ua.match(/Windows NT ([\d.]+)/))) os = 'Windows ' + ({ '10.0': '10/11', '6.3': '8.1', '6.2': '8', '6.1': '7' }[m[1]] || m[1]);
    else if ((m = ua.match(/Android ([\d.]+)/))) os = 'Android ' + m[1];
    else if ((m = ua.match(/(?:iPhone|CPU) OS ([\d_]+)/))) os = (/iPad/.test(ua) ? 'iPadOS ' : 'iOS ') + m[1].replace(/_/g, '.');
    else if (/Macintosh/.test(ua)) os = tablet ? 'iPadOS' : 'macOS';
    else if (/CrOS/.test(ua)) os = 'ChromeOS';
    else if (/Linux/.test(ua)) os = 'Linux';
    var br = 'غير معروف';
    if ((m = ua.match(/Edg(?:e|A|iOS)?\/(\d+)/))) br = 'Edge ' + m[1];
    else if ((m = ua.match(/OPR\/(\d+)/))) br = 'Opera ' + m[1];
    else if ((m = ua.match(/SamsungBrowser\/(\d+)/))) br = 'Samsung Internet ' + m[1];
    else if (/FBAN|FBAV/.test(ua)) br = 'Facebook';
    else if (/Instagram/.test(ua)) br = 'Instagram';
    else if ((m = ua.match(/(?:Chrome|CriOS)\/(\d+)/))) br = 'Chrome ' + m[1];
    else if ((m = ua.match(/(?:Firefox|FxiOS)\/(\d+)/))) br = 'Firefox ' + m[1];
    else if ((m = ua.match(/Version\/(\d+).*Safari/))) br = 'Safari ' + m[1];
    var model = '';
    if ((m = ua.match(/Android [\d.]+; ([^;)]+?)(?: Build|\))/))) model = m[1].trim() === 'K' ? '' : m[1].trim();
    else if (/iPhone/.test(ua)) model = 'iPhone'; else if (/iPad/.test(ua)) model = 'iPad';
    return { dev: dev, os: os, br: br, model: model };
  }

  function ip() {
    return Promise.race([
      fetch('https://api.ipify.org?format=json', { cache: 'no-store' }).then(function (r) { return r.json(); }).then(function (j) { return j.ip || ''; }),
      new Promise(function (r) { setTimeout(function () { r(''); }, 4000); })
    ]).catch(function () { return ''; });
  }

  function whenFirebase(cb, n) {
    n = n || 0;
    if (window.firebase && firebase.apps && firebase.apps.length && firebase.firestore && firebase.auth) return cb();
    if (n < 40) setTimeout(function () { whenFirebase(cb, n + 1); }, 250);
  }

  whenFirebase(function () {
    var done = false;
    var un = firebase.auth().onAuthStateChanged(function (u) {
      // only signed-in visitors (accounts + named guests) are logged; the database never accepts anonymous writes
      if (done || !u) return; done = true; try { un(); } catch (e) { }
      var ua = navigator.userAgent || '', d = device(ua), tz = '';
      try { tz = Intl.DateTimeFormat().resolvedOptions().timeZone || ''; } catch (e) { }
      ip().then(function (addr) {
        var rec = {
          t: firebase.firestore.FieldValue.serverTimestamp(),
          app: APP,
          path: cut(location.pathname + (location.hash || ''), 200),
          ip: cut(addr, 45),
          dev: d.dev, os: cut(d.os, 40), br: cut(d.br, 40), model: cut(d.model, 60),
          scr: cut((screen.width || 0) + '×' + (screen.height || 0), 20),
          lang: cut(navigator.language || '', 20),
          tz: cut(tz, 60),
          ref: cut(document.referrer && document.referrer.indexOf(location.origin) !== 0 ? document.referrer : '', 200),
          ua: cut(ua, 400),
          uid: u ? u.uid : '',
          name: cut(u ? (u.displayName || (window.ACCT && ACCT.fs && ACCT.fs.name && ACCT.fs.name()) || '') : '', 60),
          email: cut(u && u.email ? u.email : '', 120)
        };
        firebase.firestore().collection('visits').add(rec).catch(function () { });
      });
    });
  });
})();
