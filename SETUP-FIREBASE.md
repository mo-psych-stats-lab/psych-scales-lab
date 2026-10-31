# تفعيل الحسابات بالبريد الإلكتروني (اختياري)

الموقعان (مختبر المقاييس ومختبر الإحصاء) يعملان كاملين بدون هذه الخطوات، ويُحفظ التقدم في المتصفح. بعد التفعيل يصبح للطالب حساب واحد للموقعين، ويُحفظ تقدمه على كل أجهزته، وتظهر لوحة المتصدرين، ويعمل «النموذج الأونلاين المباشر» لجمع بيانات المقاييس.

استخدم **مشروع Firebase واحدًا** للموقعين.

## 1) إنشاء المشروع
1. ادخل إلى https://console.firebase.google.com واضغط **Add project**، وسمّه مثلًا `psych-labs` (Google Analytics غير مطلوب).
2. من **Build → Authentication → Get started**، فعّل **Email/Password** (اترك Email link معطلًا).
3. من **Authentication → Settings → Password policy**: فعّل الحد الأدنى 8 أحرف، واشترط حروفًا وأرقامًا.
4. من **Authentication → Settings → Authorized domains** أضف: `mo-psych-stats-lab.github.io`.

## 2) قاعدة البيانات
1. **Build → Firestore Database → Create database** واختر **Production mode** وأقرب منطقة (مثل `europe-west`).
2. من تبويب **Rules** الصق محتوى ملف `firestore.rules` (الموجود في المستودعين، وهو نفسه) واضغط **Publish**.
   أو من سطر الأوامر: `firebase deploy --only firestore:rules`.
3. أنشئ فهرسين للوحة المتصدرين إن طلبهما Firestore عند أول تحميل (يظهر رابط مباشر في رسالة الخطأ في المتصفح): `leaderboard` مرتب تنازليًا حسب `week`، وآخر حسب `total`.

## 3) ربط الموقع
1. **Project settings → General → Your apps → Add app → Web**، وسمّه `psych-labs-web` (لا تفعّل Hosting).
2. انسخ كائن `firebaseConfig` والصقه في ملف `firebase-config.js` في **كل مستودع** بدل `null`:
   ```js
   window.FIREBASE_CONFIG = { apiKey: "…", authDomain: "…", projectId: "…", appId: "…" };
   ```
3. ارفع التعديل (commit + push)؛ سيظهر التسجيل في صفحة «حسابي» خلال دقائق.

## 4) حماية مفتاح الـ API (مهم)
مفتاح Firebase للويب ليس سرًا (يصل لكل متصفح)، لكن قيّده:
1. افتح https://console.cloud.google.com/apis/credentials في المشروع نفسه.
2. افتح المفتاح **Browser key** واختر **Application restrictions → Websites**، وأضف:
   `https://mo-psych-stats-lab.github.io/*`
3. في **API restrictions** اختر **Restrict key** وفعّل فقط: Identity Toolkit API، Token Service API، Cloud Firestore API، Firebase Installations API، Firebase App Check API.

## 5) قوالب البريد
من **Authentication → Templates** عدّل لغة الرسائل إلى العربية، واكتب اسم المرسل «مختبرات علم النفس».

## 6) (موصى به) App Check ضد الطلبات الآلية
1. أنشئ مفتاح reCAPTCHA v3 من https://www.google.com/recaptcha/admin لنطاق `mo-psych-stats-lab.github.io`.
2. في Firebase: **Build → App Check → Apps** سجّل تطبيق الويب بمفتاح reCAPTCHA v3 السري.
3. ضع مفتاح الموقع (Site key) في `firebase-config.js`: `window.FIREBASE_APPCHECK_SITE_KEY = "…";`
4. بعد التأكد من عمل الموقع، فعّل **Enforce** لـ Firestore وAuthentication.

## 7) الحد من إساءة الاستخدام والتكلفة
- الخطة المجانية (Spark) تكفي لمئات الطلاب؛ لا تضف بطاقة دفع إلا عند الحاجة.
- لا توجد صفحات أدمن في الموقع؛ وإدارة المستخدمين تتم من لوحة Firebase فقط (Authentication → Users).
- القواعد تمنع أي حساب من قراءة بيانات غيره، وتحد حجم كل مستند، وتمنع تعديل الاستجابات بعد إرسالها، وتضع حدًا زمنيًا بين عمليات الكتابة.
