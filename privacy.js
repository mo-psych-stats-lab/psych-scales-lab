/* ===== Privacy agreement (shared by both labs) =====
   Every student who takes «اعرف نفسك» or «مكتبة الاختبارات» sees this once and must accept it. Students with an
   account accept it after signing in; students without one type their full name and the site opens a guest
   session (Firebase anonymous sign-in) so their results reach the lecturer too. Written consent is collected
   on paper. Nothing is published or shared with anyone else. ===== */
(function () {
  'use strict';
  var KEY = 'privacy-v1';
  var get = function () { try { return JSON.parse(localStorage.getItem(KEY) || '{}') || {}; } catch (e) { return {}; } };
  var setOK = function (u) { var m = get(); m[u] = Date.now(); try { localStorage.setItem(KEY, JSON.stringify(m)); } catch (e) { } };
  var fs = function () { return window.ACCT && ACCT.fs; };
  var fsOK = function () { try { return !!(fs() && ACCT.fs.ready()); } catch (e) { return false; } };
  var uid = function () { return fsOK() ? ACCT.fs.uid() : null; };
  var testsPage = function () { return /^#(know|tests)/.test(location.hash || ''); };
  function ok() { var u = uid(); return !!(u && get()[u]); }
  var POLICY = '<h3 style="margin:0" id="pvT">سياسة الخصوصية والموافقة</h3>' +
    '<p style="margin:0">قبل ما تكمل، اقرأ إزاي بنتعامل مع بياناتك:</p><ul style="margin:0;line-height:1.9;padding-inline-start:20px">' +
    '<li><b>اللي بيتحفظ:</b> اسمك، وإيميلك لو عندك حساب، وتقدمك في الدروس والتدريب، ونتايج «اعرف نفسك» و«مكتبة الاختبارات».</li>' +
    '<li><b>مين يشوفها:</b> مسؤول الموقع (المحاضر) بس، علشان المتابعة الأكاديمية وتوجيهك لو احتجت مساعدة. ومفيش طالب يقدر يشوف بيانات طالب تاني.</li>' +
    '<li><b>مش بتتنشر:</b> بياناتك مش بتتنشر ولا بتتباع ولا بتتشارك مع أي جهة تانية. ولو اتستخدمت في بحث بتكون مجمّعة ومن غير أسماء.</li>' +
    '<li><b>النتايج النفسية فرز تعليمي مش تشخيص،</b> ولو نتيجتك محتاجة اهتمام هتلاقي أرقام خطوط الدعم في صفحة النتائج.</li>' +
    '<li><b>حقك:</b> تقدر تطلب من المحاضر حذف بياناتك في أي وقت.</li></ul>';
  var BTN = 'border:0;border-radius:10px;padding:8px 18px;font:inherit;font-weight:700;background:var(--primary,var(--teal,#4f46e5));color:#fff;cursor:pointer';
  var BTN2 = 'border:1px solid currentColor;border-radius:10px;padding:8px 18px;font:inherit;background:transparent;color:inherit;cursor:pointer';
  var INP = 'border:1px solid currentColor;border-radius:9px;padding:8px 10px;font:inherit;background:transparent;color:inherit;width:100%';
  var shown = false;
  function show(guest) {
    if (shown || document.getElementById('pvOv')) return; shown = true;
    var ov = document.createElement('div'); ov.id = 'pvOv'; ov.setAttribute('role', 'dialog'); ov.setAttribute('aria-modal', 'true'); ov.setAttribute('aria-labelledby', 'pvT');
    ov.style.cssText = 'position:fixed;inset:0;z-index:9999;background:rgba(0,0,0,.45);display:grid;place-items:center;padding:16px';
    ov.innerHTML = '<div dir="rtl" style="max-width:560px;width:100%;max-height:90vh;overflow:auto;background:var(--card,var(--surface,#fff));color:var(--ink,inherit);border-radius:16px;padding:18px;display:flex;flex-direction:column;gap:10px;line-height:1.8;box-shadow:0 10px 40px rgba(0,0,0,.25)">' + POLICY +
      (guest ? '<label for="pvName" style="font-weight:700">اسمك الحقيقي بالكامل (علشان المحاضر يعرف نتيجتك)</label><input id="pvName" type="text" maxlength="40" autocomplete="name" style="' + INP + '"><p style="margin:0;font-size:13px;opacity:.8">مش لازم تعمل حساب. ولو عندك حساب، <a href="./#account">سجّل دخول</a> أحسن علشان نتايجك تتحفظ على كل أجهزتك.</p>' : '') +
      '<label style="display:flex;gap:8px;align-items:flex-start"><input type="checkbox" id="pvChk"> <span>قريت وموافق على سياسة الخصوصية دي.</span></label>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap"><button type="button" id="pvYes" style="' + BTN + '">موافق</button><button type="button" id="pvNo" style="' + BTN2 + '">مش دلوقتي</button></div><p id="pvMsg" style="margin:0;font-size:13px"></p></div>';
    document.body.appendChild(ov);
    var msg = ov.querySelector('#pvMsg'), first = ov.querySelector(guest ? '#pvName' : '#pvChk'); first.focus();
    ov.querySelector('#pvYes').onclick = async function () {
      if (guest) { var nm = ov.querySelector('#pvName').value.trim(); if (nm.split(/\s+/).length < 2) { msg.textContent = 'اكتب اسمك الأول واسم الأب على الأقل.'; return; } }
      if (!ov.querySelector('#pvChk').checked) { msg.textContent = 'علّم على «قريت وموافق» الأول.'; return; }
      if (guest) {
        var b = this; b.disabled = true; msg.textContent = 'جارٍ الحفظ…';
        try { var u = await ACCT.fs.guest(nm); setOK(u); }
        catch (e) { b.disabled = false; var c = (e && (e.code || e.message)) || ''; msg.textContent = /operation-not-allowed|admin-restricted/.test(c) ? 'الدخول من غير حساب مش متفعّل لسه في الموقع. بلّغ المحاضر، أو اعمل حساب من «حسابي».' : (e && e.message) || 'حصلت مشكلة، جرّب تاني.'; return; }
      } else setOK(uid());
      var uu = uid(); if (uu) fired[uu] = 1;
      ov.remove(); setTimeout(function () { window.dispatchEvent(new Event('privacy-accepted')); }, 800);
    };
    ov.querySelector('#pvNo').onclick = function () { ov.remove(); };
  }
  var last = null, fired = {}, guestAsked = false;
  (function poll(n) {
    var u = uid();
    if (u !== last) { last = u; shown = false; }
    if (u && !location.hash.startsWith('#admin')) { if (!get()[u]) show(false); else if (!fired[u]) { fired[u] = 1; window.dispatchEvent(new Event('privacy-accepted')); } }
    else if (!u && testsPage() && !guestAsked && fs() && ACCT.fs.canGuest && ACCT.fs.canGuest()) { guestAsked = true; shown = false; show(true); }
    if (!testsPage()) guestAsked = false;
    setTimeout(function () { poll(n + 1); }, n < 20 ? 1500 : 4000);
  })(0);
  document.addEventListener('click', function (e) { var t = e.target.closest && e.target.closest('[data-pvshow]'); if (t) { e.preventDefault(); shown = false; show(!uid() && !!(fs() && ACCT.fs.canGuest && ACCT.fs.canGuest())); } });
  window.PRIVACY = { ok: ok, show: function () { shown = false; show(!uid()); } };
})();
