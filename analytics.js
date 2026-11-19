/* ===== Admin analytics (shared by both labs; loaded only inside the admin page) =====
   Study progress of registered students (users/{uid}, readable by admins under the Firestore rules),
   daily activity, per-lesson completion, a sortable table with CSV export, and anonymous visitor totals
   from GoatCounter. Psychological results (know-yourself, tests) are not read here. ===== */
(function () {
  'use strict';
  var E = function (s) { return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); };
  var parse = function (s) { try { return JSON.parse(s || 'null'); } catch (e) { return null; } };
  var dkey = function (d) { var x = new Date(d); return x.getFullYear() + '-' + String(x.getMonth() + 1).padStart(2, '0') + '-' + String(x.getDate()).padStart(2, '0'); };
  var DAY = 864e5;
  var CSS = '.an{display:flex;flex-direction:column;gap:14px}.an-k{display:grid;grid-template-columns:repeat(auto-fill,minmax(140px,1fr));gap:8px}.an-k div{border:1px solid var(--line,#ddd);border-radius:12px;padding:10px}.an-k small{display:block;opacity:.7;font-size:12.5px}.an-k b{font-size:22px}' +
    '.an svg{width:100%;height:auto;display:block}.an-bars{display:flex;flex-direction:column;gap:4px;font-size:13px}.an-bars div{display:grid;grid-template-columns:minmax(120px,38%) 1fr 44px;gap:8px;align-items:center}.an-bars i{display:block;height:10px;border-radius:5px;background:var(--primary,var(--teal,#4f46e5))}' +
    '.an table{border-collapse:collapse;width:100%;font-size:13px}.an th,.an td{border-bottom:1px solid var(--line,#ddd);padding:6px;text-align:start;white-space:nowrap}.an th{cursor:pointer;user-select:none}.an-tw{overflow-x:auto}.an .muted{opacity:.7}';
  function css() { if (document.getElementById('an-css')) return; var s = document.createElement('style'); s.id = 'an-css'; s.textContent = CSS; document.head.appendChild(s); }
  function app() { return document.title.indexOf('الإحصاء') >= 0 ? 'stats' : 'scales'; }
  function lessonsMeta() {
    try { if (app() === 'stats' && typeof LESSONS !== 'undefined') return LESSONS.map(function (l) { return { id: l.id, t: l.t }; }); } catch (e) { }
    try { if (typeof C !== 'undefined' && C.parts) return Object.keys(C.parts).map(function (k) { return { id: k, t: (C.parts[k].title || k).split(' — ')[0], n: (C.parts[k].sections || []).length }; }); } catch (e) { }
    return [];
  }
  function row(u, A) {
    var s = parse(u[A] && u[A].state) || {}, upd = u.updated && u.updated.toDate ? u.updated.toDate() : null, days = u.days || {}, now = Date.now();
    var act30 = Object.keys(days).filter(function (k) { return now - new Date(k).getTime() <= 30 * DAY; }).length;
    var r = { id: u.id, name: u.name || '—', email: u.email || '', last: upd, xp: s.xp | 0, streak: s.streak | 0, acc: s.answered ? Math.round(100 * (s.correct || 0) / s.answered) : null, answered: s.answered | 0, labs: Object.keys(s.labs || {}).length, act30: act30, days: days, s: s };
    if (A === 'stats') { r.prog = Object.keys(s.done || {}).length; r.exams = (s.exams || []).length; r.best = (s.exams || []).reduce(function (m, e) { return Math.max(m, e.n ? Math.round(100 * e.s / e.n) : 0); }, 0); }
    else { r.prog = Object.keys(s.read || {}).length; var q = s.quiz || {}, ks = Object.keys(q); r.exams = ks.length; r.best = ks.length ? Math.round(100 * ks.reduce(function (m, k) { return m + ((q[k] && q[k].best) || 0); }, 0) / ks.length) : 0; }
    return r;
  }
  function barsDaily(R) {
    var W = 600, H = 150, n = 30, now = Date.now(), cnt = [];
    for (var i = n - 1; i >= 0; i--) { var k = dkey(now - i * DAY); cnt.push([k, R.filter(function (r) { return r.days[k]; }).length]); }
    var mx = Math.max(1, Math.max.apply(null, cnt.map(function (c) { return c[1]; }))), bw = (W - 20) / n;
    return '<svg viewBox="0 0 ' + W + ' ' + H + '" direction="ltr" role="img" aria-label="الطلاب النشطون يوميًا">' + cnt.map(function (c, i) { var h = 110 * c[1] / mx; return '<rect x="' + (10 + i * bw + 1) + '" y="' + (120 - h) + '" width="' + (bw - 2) + '" height="' + h + '" rx="2" fill="currentColor" opacity=".55"><title>' + c[0] + ': ' + c[1] + '</title></rect>' + (c[1] ? '<text x="' + (10 + i * bw + bw / 2) + '" y="' + (116 - h) + '" font-size="10" text-anchor="middle" fill="currentColor">' + c[1] + '</text>' : ''); }).join('') + '<text x="10" y="142" font-size="11" fill="currentColor" opacity=".7">' + cnt[0][0] + '</text><text x="' + (W - 10) + '" y="142" font-size="11" text-anchor="end" fill="currentColor" opacity=".7">اليوم</text></svg>';
  }
  function lessonBars(R, A) {
    var L = lessonsMeta(); if (!L.length || !R.length) return '<p class="muted">لا توجد بيانات بعد.</p>';
    var rows = L.map(function (l) {
      var c = A === 'stats' ? R.filter(function (r) { return r.s.done && r.s.done[l.id]; }).length
        : R.filter(function (r) { var rd = r.s.read || {}; return Object.keys(rd).some(function (k) { return k.indexOf(l.id + '.') === 0; }); }).length;
      return [l.t, c];
    });
    return '<div class="an-bars">' + rows.map(function (x) { var p = Math.round(100 * x[1] / R.length); return '<div><span>' + E(x[0]) + '</span><span><i style="width:' + p + '%"></i></span><span>' + p + '%</span></div>'; }).join('') + '</div>';
  }
  async function visitors() {
    var base = 'https://mo-psych-stats-lab.goatcounter.com', site = app() === 'stats' ? '/psych-stats-lab/' : '/psych-scales-lab/';
    async function get(p) { var r = await fetch(base + '/counter/' + encodeURIComponent(p) + '.json', { referrerPolicy: 'no-referrer' }); if (!r.ok) throw 0; return (await r.json()).count; }
    try { var t = await get('TOTAL'), s = await get(site).catch(function () { return '—'; }); return '<p>الزيارات منذ بداية العد: <b>' + E(t) + '</b> للموقعين معًا، منها <b>' + E(s) + '</b> للصفحة الرئيسية لهذا الموقع.</p>'; }
    catch (e) { return '<p class="muted">تعذر قراءة العدد من GoatCounter. لعرضه هنا: ادخل لوحة GoatCounter ← Settings ← فعّل «Allow adding visitor counts on your website». أو افتح اللوحة نفسها من الرابط تحت.</p>'; }
  }
  function csv(R, A) {
    var head = ['الاسم', 'الإيميل', 'آخر نشاط', A === 'stats' ? 'الدروس المكتملة' : 'الأقسام المقروءة', 'الدقة %', 'الأسئلة المجابة', 'النقاط', A === 'stats' ? 'الامتحانات' : 'اختبارات الوحدات', 'أفضل نتيجة %', 'المعامل المستخدمة', 'أيام النشاط (30 يومًا)'];
    var q = function (v) { v = String(v == null ? '' : v); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
    var lines = [head].concat(R.map(function (r) { return [r.name, r.email, r.last ? r.last.toISOString().slice(0, 16).replace('T', ' ') : '', r.prog, r.acc == null ? '' : r.acc, r.answered, r.xp, r.exams, r.best, r.labs, r.act30]; }));
    var blob = new Blob(['﻿' + lines.map(function (l) { return l.map(q).join(','); }).join('\n')], { type: 'text/csv;charset=utf-8' });
    var a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'students-' + A + '-' + dkey(Date.now()) + '.csv'; a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); }, 2000);
  }
  async function render(el, fs) {
    css(); var A = app();
    el.innerHTML = '<section class="ad-card an"><h3>تحليلات الطلاب والزوار</h3><p class="muted">جارٍ التحميل…</p></section>';
    var users = [];
    try { users = await fs.list('users'); }
    catch (e) { el.firstChild.innerHTML = '<h3>تحليلات الطلاب والزوار</h3><p style="color:#b3261e">تعذر قراءة بيانات الطلاب. انشر ملف القواعد الجديد في Firebase (القسم الخاص بـ users يسمح للمسؤول بالقراءة)، ثم حدّث الصفحة.</p>'; return; }
    var R = users.filter(function (u) { return u[A] && u[A].state; }).map(function (u) { return row(u, A); }), now = Date.now();
    var within = function (d) { return R.filter(function (r) { return r.last && now - r.last.getTime() <= d * DAY; }).length; };
    var avg = function (f) { var v = R.map(f).filter(function (x) { return x != null; }); return v.length ? Math.round(v.reduce(function (a, b) { return a + b; }, 0) / v.length) : 0; };
    var sortK = 'last', dir = -1;
    var K = function (l, v) { return '<div><small>' + l + '</small><b>' + v + '</b></div>'; };
    var sec = el.firstChild;
    sec.innerHTML = '<h3>تحليلات الطلاب والزوار</h3>' +
      '<p class="muted" style="margin:0">التقدم الدراسي لكل الطلاب المسجلين (الدروس والتدريب والنقاط وأيام النشاط). نتائج «اعرف نفسك» والاختبارات النفسية تظهر في جدولها تحت لمن وافق على مشاركتها فقط.</p>' +
      '<div class="an-k">' + K('طلاب مسجلون', R.length) + K('نشطون اليوم', within(1)) + K('آخر 7 أيام', within(7)) + K('آخر 30 يومًا', within(30)) + K(A === 'stats' ? 'متوسط الدروس المكتملة' : 'متوسط الأقسام المقروءة', avg(function (r) { return r.prog; })) + K('متوسط الدقة', avg(function (r) { return r.acc; }) + '%') + K('مجموع الأسئلة المجابة', R.reduce(function (a, r) { return a + r.answered; }, 0)) + '</div>' +
      '<h4 style="margin:6px 0 0">الطلاب النشطون يوميًا (آخر 30 يومًا)</h4>' + barsDaily(R) +
      '<h4 style="margin:6px 0 0">' + (A === 'stats' ? 'نسبة من أكملوا كل درس' : 'نسبة من بدأوا كل وحدة') + '</h4>' + lessonBars(R, A) +
      '<div style="display:flex;justify-content:space-between;align-items:center;gap:8px;flex-wrap:wrap"><h4 style="margin:6px 0 0">كل الطلاب</h4><button type="button" class="ad-btn" id="anCsv">تنزيل Excel (CSV)</button></div><div class="an-tw" id="anT"></div>' +
      '<h4 style="margin:6px 0 0">زوار الموقع (بدون تسجيل)</h4><div id="anV"><p class="muted">جارٍ القراءة…</p></div><p style="margin:0"><a href="https://mo-psych-stats-lab.goatcounter.com" target="_blank" rel="noopener">افتح لوحة GoatCounter الكاملة</a>: الصفحات والدول والأجهزة والمتصفحات ومصادر الزيارات يومًا بيوم.</p>';
    var COLS = [['name', 'الاسم'], ['email', 'الإيميل'], ['last', 'آخر نشاط'], ['prog', A === 'stats' ? 'دروس' : 'أقسام'], ['acc', 'الدقة'], ['answered', 'أسئلة'], ['xp', 'نقاط'], ['exams', A === 'stats' ? 'امتحانات' : 'اختبارات'], ['best', 'أفضل %'], ['labs', 'معامل'], ['act30', 'أيام نشاط']];
    function table() {
      var S = R.slice().sort(function (a, b) { var x = a[sortK], y = b[sortK]; if (x instanceof Date) x = x.getTime(); if (y instanceof Date) y = y.getTime(); if (x == null) x = -1; if (y == null) y = -1; return (x > y ? 1 : x < y ? -1 : 0) * dir; });
      sec.querySelector('#anT').innerHTML = R.length ? '<table><thead><tr>' + COLS.map(function (c) { return '<th data-k="' + c[0] + '">' + c[1] + (sortK === c[0] ? (dir < 0 ? ' ▼' : ' ▲') : '') + '</th>'; }).join('') + '</tr></thead><tbody>' + S.map(function (r) { return '<tr><td><b>' + E(r.name) + '</b></td><td dir="ltr">' + (r.email ? E(r.email) : '<span dir="rtl" class="muted">ضيف</span>') + '</td><td>' + (r.last ? r.last.toLocaleString('ar-EG', { dateStyle: 'short', timeStyle: 'short' }) : '—') + '</td><td>' + r.prog + '</td><td>' + (r.acc == null ? '—' : r.acc + '%') + '</td><td>' + r.answered + '</td><td>' + r.xp + '</td><td>' + r.exams + '</td><td>' + (r.exams ? r.best + '%' : '—') + '</td><td>' + r.labs + '</td><td>' + r.act30 + '</td></tr>'; }).join('') + '</tbody></table>' : '<p class="muted">لا يوجد طلاب مسجلون بعد.</p>';
      sec.querySelectorAll('#anT th').forEach(function (th) { th.onclick = function () { if (sortK === th.dataset.k) dir = -dir; else { sortK = th.dataset.k; dir = -1; } table(); }; });
    }
    table();
    sec.querySelector('#anCsv').onclick = function () { csv(R, A); };
    sec.querySelector('#anV').innerHTML = await visitors();
  }
  window.ADMIN_AN = { render: render };
})();
